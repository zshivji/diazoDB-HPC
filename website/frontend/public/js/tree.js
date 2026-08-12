console.log("tree.js is running")

// define variables
const tooltip = document.getElementById("tooltip")
const treeContainer = document.getElementById("tree-container")
const searchInput = document.getElementById("tree-search")
const searchExampleButtons = document.querySelectorAll("[data-search-example]")
let labelData = []
let metadataKeyLookup = new Map()

// define colors for each gene
const GENE_COLORS = {
  nifA: "#c084fc", // purple
  nifB: "#f472b6", // pink
  nifH: "#60a5fa", // blue
  nifD: "#f59e0b", // orange
  nifK: "#10b981", // green
  nifE: "#f97316", // amber
  nifN: "#14b8a6", // teal
  nifV: "#a3e635", // lime
  anfG: "#818cf8", // indigo
  anfO: "#fb7185", // rose
  vnfG: "#22d3ee", // cyan
}

const DEFAULT_GENE_COLOR = "#9ca3af"
const SHOW_OPERON_IN_TOOLTIP = true

// helper functions

function escapeHtml(str) {
  return String(str == null ? "" : str)
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;")
    .replace(/'/g, "&#039;")
}

function clamp(value, min, max) {
  return Math.max(min, Math.min(max, value))
}

function formatBp(bp) {
  if (bp >= 1000000) return `${(bp / 1000000).toFixed(1)} Mb`
  if (bp >= 1000) return `${(bp / 1000).toFixed(1)} kb`
  return `${Math.round(bp)} bp`
}

function normalizeDirection(direction) {
  if (Number(direction) < 0) return "reverse"

  const d = String(direction || "").toLowerCase()
  if (d === "reverse" || d === "-" || d === "-1" || d === "left") {
    return "reverse"
  }
  return "forward"
}

function getGeneColor(geneName) {
  const key = String(geneName || "").trim()
  return GENE_COLORS[key] || DEFAULT_GENE_COLOR
}

function normalizeTreeIdentifier(value) {
  const normalized = String(value || "")
    .trim()
    .replace(/^'+|'+$/g, "")
    .replace(/\b(GB|RS)\s+(GCA|GCF)\s+(\d+\.\d+)\b/gi, "$1_$2_$3")
    .replace(/\b([A-Z]{2})\s+([A-Z]*\d+(?:\.\d+)?)\b/g, "$1_$2")

  // Tree labels and metadata keys use the same fields but may differ in
  // whitespace around separators (and in accession formatting). Keep the
  // final copy number: it is part of the exact metadata ID.
  return normalized
    .split("|")
    .map((part) => part.trim())
    .filter(Boolean)
    .join(" | ")
}

function canonicalLookupKey(value) {
  const parts = normalizeTreeIdentifier(value).split("|")
  if (parts.length > 0) {
    parts[0] = parts[0].replace(/_/g, " ")
  }
  return parts.join("|")
}

function isGenomeId(value) {
  return /^(GB|RS)_(GCA|GCF)_\d+\.\d+$/i.test(String(value || "").trim())
}

function getGenomeContigLookupKey(value) {
  const parts = normalizeTreeIdentifier(value).split("|")
  const genomeIndex = parts.findIndex(isGenomeId)

  if (genomeIndex === -1 || genomeIndex + 1 >= parts.length) {
    return ""
  }

  return `${parts[genomeIndex]}|${parts[genomeIndex + 1]}`
}

function addLookupValue(lookup, key, value) {
  if (key && !lookup.has(key)) {
    lookup.set(key, value)
  }
}

function createMetadataLookup(metadata) {
  const lookup = new Map()

  Object.keys(metadata).forEach((key) => {
    addLookupValue(lookup, normalizeTreeIdentifier(key), key)
    addLookupValue(lookup, canonicalLookupKey(key), key)
    addLookupValue(lookup, getGenomeContigLookupKey(key), key)
  })

  return lookup
}

function findMetadataKey(candidate, metadata) {
  if (!candidate) return null
  if (metadata[candidate]) return candidate

  const lookupKeys = [
    normalizeTreeIdentifier(candidate),
    canonicalLookupKey(candidate),
    getGenomeContigLookupKey(candidate),
  ]

  for (const lookupKey of lookupKeys) {
    const metadataKey = metadataKeyLookup.get(lookupKey)
    if (metadataKey && metadata[metadataKey]) return metadataKey
  }

  return null
}

// render operon

function renderOperonHTML(nodeMeta) {
  if (!nodeMeta || !nodeMeta.operon || !nodeMeta.operon.genes || nodeMeta.operon.genes.length === 0) {
    return ""
  }

  const genes = [...nodeMeta.operon.genes].sort(
    (a, b) => Number(a.start) - Number(b.start),
  )

  let regionStart = Number(nodeMeta.operon.region_start)
  let regionEnd = Number(nodeMeta.operon.region_end)

  if (Number.isNaN(regionStart)) {
    regionStart = Math.min(...genes.map((g) => Number(g.start)))
  }
  if (Number.isNaN(regionEnd)) {
    regionEnd = Math.max(...genes.map((g) => Number(g.end)))
  }

  const span = Math.max(1, regionEnd - regionStart)
  const leftPad = 8
  const rightPad = 8
  const trackWidth = 404 // approx tooltip inner width minus padding
  const usableWidth = trackWidth - leftPad - rightPad

  const geneHtml = genes
    .map((gene) => {
      const start = Number(gene.start)
      const end = Number(gene.end)
      const direction = normalizeDirection(gene.direction)
      const color = getGeneColor(gene.gene_name)

      const leftPx = leftPad + ((start - regionStart) / span) * usableWidth
      const rightPx = leftPad + ((end - regionStart) / span) * usableWidth
      const widthPx = Math.max(10, rightPx - leftPx)
      const geneName = String(gene.gene_name || "").trim()
      const label = escapeHtml(geneName || gene.gene_id || "")
      const displayLabel = /^WP_/.test(geneName) ? "" : label
      const product = escapeHtml(gene.product || "")
      const title = `${label} | ${direction} | ${start}-${end}${product ? ` | ${product}` : ""}`

      if (direction === "reverse") {
        return `
        <div class="tooltip-operon-gene"
             style="left:${leftPx}px; width:${widthPx}px;"
             title="${title}">
          <div class="tooltip-operon-gene-label">${displayLabel}</div>
          <div class="tooltip-operon-arrow-reverse"
               style="border-right-color:${color};"></div>
          <div class="tooltip-operon-gene-body"
               style="background-color:${color};"></div>
        </div>
      `
      }

      return `
      <div class="tooltip-operon-gene"
           style="left:${leftPx}px; width:${widthPx}px;"
           title="${title}">
        <div class="tooltip-operon-gene-label">${displayLabel}</div>
        <div class="tooltip-operon-gene-body"
             style="background-color:${color};"></div>
        <div class="tooltip-operon-arrow-forward"
             style="border-left-color:${color};"></div>
      </div>
    `
    })
    .join("")

  const scaleHtml = renderScaleHTML(regionStart, regionEnd)

  return `
    <div class="tooltip-operon-wrapper">
      <div class="tooltip-operon-track">
        ${geneHtml}
      </div>
      ${scaleHtml}
    </div>
  `
}

function renderScaleHTML(regionStart, regionEnd) {
  const fractions = [0, 0.2, 0.4, 0.6, 0.8, 1]

  const ticks = fractions
    .map((fraction, i) => {
      const value = Math.round(
        regionStart + fraction * (regionEnd - regionStart),
      )
      const label = formatBp(value)

      let extraClass = ""
      if (i === 0) extraClass = "start"
      if (i === fractions.length - 1) extraClass = "end"

      return `
      <div class="tooltip-operon-scale-tick ${extraClass}" style="left:${fraction * 100}%;">
        <div class="tooltip-operon-scale-mark"></div>
        <span class="tooltip-operon-scale-label">${label}</span>
      </div>
    `
    })
    .join("")

  return `
    <div class="tooltip-operon-scale">
      <div class="tooltip-operon-scale-line"></div>
      ${ticks}
    </div>
  `
}

function resolveMetadataKey(rawLabel, metadata) {
  const rawId = normalizeTreeIdentifier(rawLabel)
  const normalizedPipeKey = normalizeTreeIdentifier(rawId)

  const withoutLeafSuffix = rawId.replace(/([_\s])\d+$/, "")
  const candidates = [rawId, withoutLeafSuffix, normalizedPipeKey]

  for (const candidate of candidates) {
    const metadataKey = findMetadataKey(candidate, metadata)
    if (metadataKey) return metadataKey
  }

  // Legacy fallback: build accession-style keys like GB_GCA_031257675.1_nif.
  const parts = normalizedPipeKey.split("|")
  for (const part of parts) {
    const match = part.trim().match(/\b(GB|RS)[ _](GCA|GCF)[ _](\d+\.\d+)\b/i)
    if (!match) continue

    const accessionKey = `${match[1].toUpperCase()}_${match[2].toUpperCase()}_${match[3]}_nif`
    const metadataKey = findMetadataKey(accessionKey, metadata)
    if (metadataKey) return metadataKey
  }

  return null
}

function normalizeText(value, fallback = "Not available") {
  const text = String(value == null ? "" : value).trim()
  return text ? text : fallback
}

function normalizeEnvironment(value) {
  const text = String(value == null ? "" : value).trim()
  return !text || text.toLowerCase() === "none" ? "Not Available" : text
}

function getTreeIdParts(id) {
  return String(id == null ? "" : id)
    .split("|")
    .map((part) => part.trim())
}

function splitTreeLabelText(value) {
  const fullValue = String(value == null ? "" : value).trim()
  if (!fullValue) {
    return { primary: "", secondary: "" }
  }

  const parts = fullValue
    .split("|")
    .map((part) => part.trim())
    .filter(Boolean)

  // The first two fields are always the organism and nitrogenase type
  // (nif/vnf/anf). Everything after them is the metadata ID suffix.
  if (parts.length >= 3) {
    return {
      primary: parts.slice(0, 2).join(" | "),
      secondary: ` | ${parts.slice(2).join(" | ")}`,
    }
  }

  return { primary: fullValue, secondary: "" }
}

function renderTreeLabelWithMetadataStyle(label, value) {
  const { primary, secondary } = splitTreeLabelText(value)
  if (!primary && !secondary) return

  const svgNs = "http://www.w3.org/2000/svg"
  label.textContent = ""

  const primarySpan = document.createElementNS(svgNs, "tspan")
  primarySpan.setAttribute("fill", "#111827")
  primarySpan.setAttribute("font-weight", "600")
  primarySpan.textContent = primary
  label.appendChild(primarySpan)

  if (secondary) {
    const secondarySpan = document.createElementNS(svgNs, "tspan")
    secondarySpan.setAttribute("fill", "#ffffff")
    secondarySpan.setAttribute("font-weight", "500")
    secondarySpan.textContent = secondary
    label.appendChild(secondarySpan)
  }
}

function formatTaxonomy(taxonomy) {
  const raw = String(taxonomy == null ? "" : taxonomy).trim()
  if (!raw) return "Not available"

  return raw
    .split(";")
    .map((part) => part.trim().replace(/^[a-z]__/, ""))
    .filter(Boolean)
    .join(" > ")
}

function toSearchText(parts) {
  return parts
    .map((value) =>
      String(value == null ? "" : value)
        .toLowerCase()
        .trim(),
    )
    .filter(Boolean)
    .join(" ")
}

function buildTooltipHTML(id, nodeMeta) {
  const idParts = getTreeIdParts(id)
  const title = idParts.length >= 2
    ? `${idParts[0]} | ${idParts[1]}`
    : normalizeText(nodeMeta.organism || nodeMeta.species || id, id)
  const genome = normalizeText(nodeMeta.genome)
  const contig = normalizeText(idParts[3])
  const group = normalizeText(nodeMeta.group)
  const environments = normalizeEnvironment(
    nodeMeta.environments || nodeMeta.environment,
  )
  const regulon = normalizeText(nodeMeta.regulon)
  const taxonomy = formatTaxonomy(nodeMeta.taxonomy)
  const operonHtml = SHOW_OPERON_IN_TOOLTIP ? renderOperonHTML(nodeMeta) : ""

  return `
    <div class="tooltip-card">
      <div class="tooltip-title">${escapeHtml(title)}</div>

      <div class="tooltip-meta-grid">
        <div class="tooltip-key">Genome</div>
        <div class="tooltip-value">${escapeHtml(genome)}</div>

        <div class="tooltip-key">Contig</div>
        <div class="tooltip-value">${escapeHtml(contig)}</div>

        <div class="tooltip-key">Group</div>
        <div class="tooltip-value">${escapeHtml(group)}</div>

        <div class="tooltip-key">Environment</div>
        <div class="tooltip-value">${escapeHtml(environments)}</div>

        <div class="tooltip-key">Regulon</div>
        <div class="tooltip-value">${escapeHtml(regulon)}</div>
      </div>

      <div class="tooltip-taxonomy-label">Taxonomy</div>
      <div class="tooltip-taxonomy">${escapeHtml(taxonomy)}</div>

      ${operonHtml}
    </div>
  `
}

function getOrganismDetailUrl(id) {
  return `/organism?id=${encodeURIComponent(id)}`
}

function fetchMetadata(url) {
  return fetch(url).then((response) => {
    if (!response.ok) {
      throw new Error(`Failed to load metadata.json: ${response.status} ${response.statusText}`)
    }

    return response.text().then((text) => {
      // The generated metadata contains Python-style NaN values, which are
      // not valid JSON and cause response.json() to fail in the browser.
      const validJson = text.replace(/\bNaN\b/g, "null")
      return JSON.parse(validJson)
    })
  })
}

function positionTooltip(e) {
  const offset = 12

  tooltip.style.left = "0px"
  tooltip.style.top = "0px"
  tooltip.style.display = "block"

  const width = tooltip.offsetWidth
  const height = tooltip.offsetHeight

  let left = e.pageX + offset
  let top = e.pageY + offset

  if (left + width > window.scrollX + window.innerWidth - 12) {
    left = e.pageX - width - offset
  }
  if (top + height > window.scrollY + window.innerHeight - 12) {
    top = e.pageY - height - offset
  }

  left = clamp(
    left,
    window.scrollX + 8,
    window.scrollX + window.innerWidth - width - 8,
  )
  top = clamp(
    top,
    window.scrollY + 8,
    window.scrollY + window.innerHeight - height - 8,
  )

  tooltip.style.left = `${left}px`
  tooltip.style.top = `${top}px`
}

// main

if (treeContainer) {
  const svgUrl = new URL("./images/nifH-tree.svg", window.location.href).href
  const metadataUrl = new URL("./results/metadata.json", window.location.href).href

  // load SVG + metadata + leaf labels + add tooltip interactivity
  Promise.all([
    fetch(svgUrl).then((r) => {
      if (!r.ok) {
        throw new Error(`Failed to load nifH-tree.svg: ${r.status} ${r.statusText}`)
      }
      return r.text()
    }),
    fetchMetadata(metadataUrl),
  ])
    .then(([svg, metadata]) => {
      treeContainer.innerHTML = svg
      metadataKeyLookup = createMetadataLookup(metadata)
      console.log("SVG + metadata loaded")

    // grab iTOL leaf labels (<text> elements)
    const labels = document.querySelectorAll("svg text")
    console.log("labels:", labels.length)

    labelData = [] // reset labelData

    labels.forEach((label) => {
      const rawId = label.textContent.trim().replace(/'/g, "") // remove quotes
      const id = resolveMetadataKey(rawId, metadata)

      if (!id || !metadata[id]) return // skip label if not in metadata

      // store label + metadata for search functionality
      const nodeMeta = metadata[id]
      labelData.push({
        label,
        id,
        text: toSearchText([
          rawId,
          id,
          nodeMeta.organism,
          nodeMeta.species,
          nodeMeta.genome,
          nodeMeta.group,
          nodeMeta.environments,
          nodeMeta.environment,
          nodeMeta.regulon,
        ]),
      })

      label.classList.add("tree-label")
      renderTreeLabelWithMetadataStyle(label, id || rawId)
      label.setAttribute("role", "link")
      label.setAttribute("tabindex", "0")
      label.setAttribute(
        "aria-label",
        `Open organism page for ${nodeMeta.organism || nodeMeta.species || id}`,
      )

      label.addEventListener("mouseover", (e) => {
        tooltip.innerHTML = buildTooltipHTML(id, metadata[id])
        tooltip.style.display = "block"
        positionTooltip(e)
      })

      label.addEventListener("mousemove", (e) => {
        positionTooltip(e)
      })

      label.addEventListener("mouseout", () => {
        tooltip.style.display = "none"
      })

      label.addEventListener("click", () => {
        window.location.href = getOrganismDetailUrl(id)
      })

      label.addEventListener("keydown", (e) => {
        if (e.key !== "Enter" && e.key !== " ") return
        e.preventDefault()
        window.location.href = getOrganismDetailUrl(id)
      })
    })

      console.log("metadata-linked labels:", labelData.length)
    })
    .catch((error) => {
      console.error("Failed to initialize tree:", error)
    })
}

// add search functionality
if (searchInput) {
  searchExampleButtons.forEach((button) => {
    button.addEventListener("click", () => {
      searchInput.value = button.dataset.searchExample
      searchInput.dispatchEvent(new Event("input", { bubbles: true }))
      searchInput.focus()
    })
  })

  searchInput.addEventListener("input", () => {
    const query = searchInput.value.toLowerCase().trim()
    const tokens = query.split(/\s+/).filter(Boolean)

    if (!query) {
      labelData.forEach((d) => {
        d.label.classList.remove("tree-match", "tree-dim")
      })
      return
    }

    labelData.forEach((d) => {
      const isMatch = tokens.every((token) => d.text.includes(token))

      if (isMatch) {
        d.label.classList.add("tree-match")
        d.label.classList.remove("tree-dim")
      } else {
        d.label.classList.remove("tree-match")
        d.label.classList.add("tree-dim")
      }
    })
  })
}
