#!/usr/bin/env bash
# Robust AlphaFold-first fetcher from <RefSeq>\t<UniProt|MISSING>
# Order: AF files (v6→v4) → AF API → AF entry HTML → AF search mapping → UniProt primary → SWISS-MODEL
# Bash-3.2 friendly (no ^^, jq, or mapfile). Quiet curl (no 404 spam).

set -euo pipefail
IFS=$'\n\t'

INPUT="${1:?Input TSV (RefSeq<TAB>UniProt|MISSING) required}"
OUT="${2:?Output directory required}"
mkdir -p "$OUT"

: "${SWM_MAX:=1}"
: "${CURL_CONNECT_TIMEOUT:=10}"
: "${CURL_MAX_TIME:=60}"
: "${CURL_RETRIES:=3}"

AF_HOST="https://alphafold.ebi.ac.uk"
AF_FILES="${AF_HOST}/files"
AF_API="${AF_HOST}/api/prediction"
AF_SEARCH="${AF_HOST}/api/search"
SWM_HOST="https://swissmodel.expasy.org"
SWM_PAGE_BASE="${SWM_HOST}/repository/uniprot"
UA="Mozilla/5.0 (Macintosh; Intel Mac OS X) AppleWebKit/537.36 (KHTML, like Gecko) curl/8 Safari/537.36"

log(){ printf '%s\n' "$*" >&2; }

# ---------- curl helpers (quiet on expected failures) ----------
curl_head_ok() {  # url
  curl -fsSI -A "$UA" -L \
    --retry "$CURL_RETRIES" --retry-all-errors --retry-delay 1 \
    --connect-timeout "$CURL_CONNECT_TIMEOUT" --max-time "$CURL_MAX_TIME" \
    "$1" >/dev/null 2>&1
}
curl_get_quiet() {  # url out
  curl -fsSL -A "$UA" -L \
    --retry "$CURL_RETRIES" --retry-all-errors --retry-delay 1 \
    --connect-timeout "$CURL_CONNECT_TIMEOUT" --max-time "$CURL_MAX_TIME" \
    -o "$2" "$1" >/dev/null 2>&1
}
curl_out_quiet() {  # url -> stdout
  curl -fsSL -A "$UA" -L \
    --retry "$CURL_RETRIES" --retry-all-errors --retry-delay 1 \
    --connect-timeout "$CURL_CONNECT_TIMEOUT" --max-time "$CURL_MAX_TIME" \
    "$1" 2>/dev/null
}

# ---------- small utils (Bash 3 safe) ----------
upper() { printf '%s' "$1" | tr '[:lower:]' '[:upper:]'; }
iso_base() { printf '%s' "$1" | sed -E 's/-[0-9]+$//'; }

# Build candidate ID list: given -> base
candidates_from_raw() {
  local raw="$1" U B
  U="$(upper "$raw")"
  B="$(iso_base "$U")"
  printf '%s\n' "$U"
  [ "$B" != "$U" ] && printf '%s\n' "$B" || true
}

# UniProt: get primary accession (and its base) for an arbitrary token
uniprot_primary_candidates() {
  local tok="$1" js tsv primary pbase
  js="$(curl_out_quiet "https://rest.uniprot.org/uniprotkb/$(upper "$tok").json" || true)"
  if [ -n "$js" ]; then
    primary="$(printf '%s' "$js" | grep -Eo '"primaryAccession"\s*:\s*"[^"]+"' | head -n1 | sed -E 's/.*"([^"]+)".*/\1/')"
  fi
  if [ -z "${primary-}" ]; then
    tsv="$(curl_out_quiet "https://rest.uniprot.org/uniprotkb/search?query=(accession:$(upper "$tok")%20OR%20id:$(upper "$tok"))&fields=primaryAccession&format=tsv&size=1" || true)"
    primary="$(printf '%s' "$tsv" | awk 'NR==2{print $1}')"
  fi
  if [ -n "${primary-}" ]; then
    printf '%s\n' "$(upper "$primary")"
    pbase="$(iso_base "$primary")"
    [ "$pbase" != "$primary" ] && printf '%s\n' "$(upper "$pbase")" || true
  fi
}

# AlphaFold search: map weird token -> canonical db_id
af_search_map() {
  local q="$1" js id
  js="$(curl_out_quiet "${AF_SEARCH}?query=$(upper "$q")" || true)"
  id="$(printf '%s' "$js" | grep -Eo '"db_id"\s*:\s*"[^"]+"' | head -n1 | sed -E 's/.*"([^"]+)".*/\1/')"
  [ -n "$id" ] && printf '%s' "$id" || return 1
}

# ---------- AlphaFold lookups ----------
# 1) Direct files probe: v6 then v4, F1..F5
af_files_try() {
  local id="$1" ver F url
  for ver in v6 v4; do
    for F in F1 F2 F3 F4 F5; do
      url="${AF_FILES}/AF-${id}-${F}-model_${ver}.pdb"
      if curl_head_ok "$url"; then
        printf '%s\n' "$url"
        return 0
      fi
    done
  done
  return 1
}

# 2) API
AF_PDB_URL=""; AF_PAE_URL=""
af_api_try() {
  local id="$1" js
  AF_PDB_URL=""; AF_PAE_URL=""
  js="$(curl_out_quiet "${AF_API}/${id}" || true)"
  [ -z "$js" ] && return 1
  AF_PDB_URL="$(printf '%s' "$js" | grep -Eo '"pdbUrl"\s*:\s*"[^"]+"' | head -n1 | sed -E 's/.*"([^"]+)".*/\1/')"
  [ -z "$AF_PDB_URL" ] && return 1
  AF_PAE_URL="$(printf '%s' "$js" | grep -Eo '"paeDownloadUrl"\s*:\s*"[^"]+"' | head -n1 | sed -E 's/.*"([^"]+)".*/\1/')"
  return 0
}

# 3) HTML scrape
AF_HTML_PDB_URL=""; AF_HTML_PAE_URL=""
af_html_try() {
  local id="$1" html
  AF_HTML_PDB_URL=""; AF_HTML_PAE_URL=""
  html="$(curl_out_quiet "${AF_HOST}/entry/${id}" || true)"
  [ -z "$html" ] && return 1
  AF_HTML_PDB_URL="$(printf '%s' "$html" | grep -Eo '/files/AF-[A-Za-z0-9]+(-[0-9]+)?-F[1-5]-model_v(6|4)\.pdb' | head -n1)"
  [ -z "$AF_HTML_PDB_URL" ] && return 1
  AF_HTML_PDB_URL="${AF_HOST}${AF_HTML_PDB_URL}"
  AF_HTML_PAE_URL="$(printf '%s' "$html" | grep -Eo '/files/AF-[A-Za-z0-9]+(-[0-9]+)?-F[1-5]-predicted_aligned_error_v[0-9]+\.json' | head -n1)"
  [ -n "$AF_HTML_PAE_URL" ] && AF_HTML_PAE_URL="${AF_HOST}${AF_HTML_PAE_URL}" || true
  return 0
}

# Download helper (logs, tries to pick PAE variant)
af_download_from_files_url() {
  local refseq="$1" files_url="$2" outp stem pae
  outp="${OUT}/${refseq}__$(basename "$files_url")"
  log "AF FILE  ${refseq}  url=$(basename "$files_url")"
  curl_get_quiet "$files_url" "$outp" || { log "AF FAIL  ${refseq}  (files)"; return 1; }
  stem="${files_url%-model_*}"
  for pae in "${stem}-predicted_aligned_error_v6.json" "${stem}-predicted_aligned_error_v4.json" "${stem}-predicted_aligned_error_v1.json"; do
    curl_head_ok "$pae" && curl_get_quiet "$pae" "${OUT}/${refseq}__$(basename "$pae")" || true
  done
  log "AF OK    ${refseq}  -> $(basename "$outp")  [files]"
  return 0
}

af_download_from_api() {
  local refseq="$1" pdb_url="$2" pae_url="${3-}" outp
  outp="${OUT}/${refseq}__$(basename "$pdb_url")"
  log "AF API   ${refseq}  url=$(basename "$pdb_url")"
  curl_get_quiet "$pdb_url" "$outp" || { log "AF FAIL  ${refseq}  (api)"; return 1; }
  [ -n "$pae_url" ] && curl_head_ok "$pae_url" && curl_get_quiet "$pae_url" "${OUT}/${refseq}__$(basename "$pae_url")" || true
  log "AF OK    ${refseq}  -> $(basename "$outp")  [api]"
  return 0
}

af_download_from_html() {
  local refseq="$1" pdb_url="$2" pae_url="${3-}" outp
  outp="${OUT}/${refseq}__$(basename "$pdb_url")"
  log "AF HTML  ${refseq}  url=$(basename "$pdb_url")"
  curl_get_quiet "$pdb_url" "$outp" || { log "AF FAIL  ${refseq}  (html)"; return 1; }
  [ -n "$pae_url" ] && curl_head_ok "$pae_url" && curl_get_quiet "$pae_url" "${OUT}/${refseq}__$(basename "$pae_url")" || true
  log "AF OK    ${refseq}  -> $(basename "$outp")  [html]"
  return 0
}

# Try AF for a single ID with the three mechanisms (files→api→html)
af_try_all_for_id() {
  local refseq="$1" id="$2" files_url
  if files_url="$(af_files_try "$id")"; then
    af_download_from_files_url "$refseq" "$files_url" && return 0
  fi
  if af_api_try "$id"; then
    af_download_from_api "$refseq" "$AF_PDB_URL" "${AF_PAE_URL-}" && return 0
  fi
  if af_html_try "$id"; then
    af_download_from_html "$refseq" "$AF_HTML_PDB_URL" "${AF_HTML_PAE_URL-}" && return 0
  fi
  return 1
}

fetch_alphafold_one() {
  local refseq="$1" raw="$2" id mapped seen=""
  # 1) direct candidates (given + isoform-stripped)
  while read -r id; do
    [ -z "$id" ] && continue
    printf '%s\n' "$seen" | grep -qx "$id" || { seen="$seen
$id"; af_try_all_for_id "$refseq" "$id" && return 0; }
  done < <(candidates_from_raw "$raw")

  # 2) AF search mapping of each candidate → canonical db_id
  while read -r id; do
    [ -z "$id" ] && continue
    mapped="$(af_search_map "$id" || true)"
    if [ -n "$mapped" ]; then
      log "AF MAP   ${refseq}  ${id} -> ${mapped}"
      printf '%s\n' "$seen" | grep -qx "$mapped" || { seen="$seen
$mapped"; af_try_all_for_id "$refseq" "$mapped" && return 0; }
    fi
  done < <(candidates_from_raw "$raw")

  # 3) UniProt primary accession and base
  while read -r id; do
    [ -z "$id" ] && continue
    printf '%s\n' "$seen" | grep -qx "$id" || { seen="$seen
$id"; af_try_all_for_id "$refseq" "$id" && return 0; }
  done < <(uniprot_primary_candidates "$raw" || true)

  log "AF NONE  ${refseq}  ${raw}"
  return 1
}

# ---------- SWISS-MODEL fallback ----------
extract_swm_pdb_links_primary() {
  grep -Eo 'https://swissmodel\.expasy\.org/repository/[A-Za-z0-9]{16,64}\.pdb' || true
  grep -Eo 'href="/repository/[A-Za-z0-9]{16,64}\.pdb"' | sed -E 's/^href="(.*)"/\1/' \
    | sed -E "s#^/repository/#${SWM_HOST}/repository/#" || true
}
extract_swm_pdb_links_fallback() {
  grep -Eo 'id="cid_[a-f0-9]{16,64}"' \
    | sed -E 's/^id="cid_([a-f0-9]{16,64})".*/\1/' \
    | awk -v host="$SWM_HOST" '{print host "/repository/" $1 ".pdb"}' || true
}
fetch_swissmodel() {
  local refseq="$1" uniprot="$2" page="${SWM_PAGE_BASE}/$(upper "$uniprot")"
  local tmp html; tmp="$(mktemp)"
  if ! curl -fsSL -A "$UA" -L \
        --retry "$CURL_RETRIES" --retry-all-errors --retry-delay 1 \
        --connect-timeout "$CURL_CONNECT_TIMEOUT" --max-time "$CURL_MAX_TIME" \
        "$page" > "$tmp" 2>/dev/null; then
    log "SWM ERR  ${refseq}  ${uniprot}  (load fail)"; rm -f "$tmp"; return 1
  fi
  html="$(cat "$tmp")"; rm -f "$tmp"
  # collect links without mapfile
  local links; links="$(printf '%s' "$html" | extract_swm_pdb_links_primary; printf '%s' "$html" | extract_swm_pdb_links_fallback)" || true
  [ -z "$links" ] && { log "SWM NA   ${refseq}  ${uniprot}  (no pdb links)"; return 1; }
  local saved=0 url
  # unique & limited
  printf '%s\n' "$links" | awk 'NF && !seen[$0]++' | while read -r url; do
    [ $saved -ge $SWM_MAX ] && break
    [ -z "$url" ] && continue
    local base coord outp
    base="$(basename "$url")"; coord="${base%.pdb}"
    outp="${OUT}/${refseq}__SWM-$(upper "$uniprot")-${coord}.pdb"
    log "SWM GET  ${refseq}  $(upper "$uniprot")  $(basename "$url")"
    if curl_get_quiet "$url" "$outp"; then
      log "SWM OK   ${refseq}  -> $(basename "$outp")"
      saved=$((saved+1))
    fi
  done
  [ "$saved" -gt 0 ]
}

# ---------- main ----------
process_line() {
  local line="$1" refseq uniprot
  [ -z "$line" ] && return 0
  refseq="$(printf '%s' "$line" | awk -F'\t' '{print $1}')"
  uniprot="$(printf '%s' "$line" | awk -F'\t' '{print $2}')"
  [ -z "$refseq" ] || [ -z "$uniprot" ] && return 0
  case "$refseq" in \#*|RefSeq*|refseq*) return 0 ;; esac
  case "$uniprot" in UniProt*|uniprot*) return 0 ;; esac

  case "$(upper "$uniprot")" in
    MISSING|NA|-) log "MISS     ${refseq}  (no UniProt ID)"; return 0 ;;
  esac

  if fetch_alphafold_one "$refseq" "$uniprot"; then
    return 0
  fi
  fetch_swissmodel "$refseq" "$uniprot" || true
}

# robust to CRLF
while IFS= read -r line || [ -n "$line" ]; do
  process_line "$(printf '%s' "$line" | tr -d '\r')"
done < "$INPUT"
