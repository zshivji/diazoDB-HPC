console.log("tree.js is running")

// define variables
const tooltip = document.getElementById("tooltip")
const searchInput = document.getElementById("tree-search")
let labelData = []

// define colors for each gene
const GENE_COLORS = {
  nifH: "#60a5fa", // blue
  nifD: "#f59e0b", // orange
  nifK: "#10b981", // green
}

const DEFAULT_GENE_COLOR = "#9ca3af"
const SHOW_OPERON_IN_TOOLTIP = false

// helper functions

function escapeHtml(str) {
  return String(str ?? "")
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
  const d = String(direction || "").toLowerCase()
  if (d === "reverse" || d === "-" || d === "left") return "reverse"
  return "forward"
}

function getGeneColor(geneName) {
  const key = String(geneName || "").trim()
  return GENE_COLORS[key] || DEFAULT_GENE_COLOR
}

// render operon

function renderOperonHTML(nodeMeta) {
  if (!nodeMeta?.operon?.genes || nodeMeta.operon.genes.length === 0) {
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
      const label = escapeHtml(gene.gene_name || gene.gene_id || "")
      const product = escapeHtml(gene.product || "")
      const title = `${label} | ${direction} | ${start}-${end}${product ? ` | ${product}` : ""}`

      if (direction === "reverse") {
        return `
        <div class="tooltip-operon-gene"
             style="left:${leftPx}px; width:${widthPx}px;"
             title="${title}">
          <div class="tooltip-operon-gene-label">${label}</div>
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
        <div class="tooltip-operon-gene-label">${label}</div>
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
  const rawId = String(rawLabel || "")
    .trim()
    .replace(/^'+|'+$/g, "")

  const normalizedPipeKey = rawId
    .replace(/\b(GB|RS)\s+(GCA|GCF)\s+(\d+\.\d+)\b/gi, "$1_$2_$3")
    .replace(/\|([^|]+?)\s+\d+$/, "|$1")

  const withoutLeafSuffix = rawId.replace(/([_\s])\d+$/, "")
  const candidates = [rawId, withoutLeafSuffix, normalizedPipeKey]

  for (const candidate of candidates) {
    if (candidate && metadata[candidate]) return candidate
  }

  // Legacy fallback: build accession-style keys like GB_GCA_031257675.1_nif.
  const parts = normalizedPipeKey.split("|")
  for (const part of parts) {
    const match = part.trim().match(/\b(GB|RS)[ _](GCA|GCF)[ _](\d+\.\d+)\b/i)
    if (!match) continue

    const accessionKey = `${match[1].toUpperCase()}_${match[2].toUpperCase()}_${match[3]}_nif`
    if (metadata[accessionKey]) return accessionKey
  }

  return null
}

function normalizeText(value, fallback = "Not available") {
  const text = String(value ?? "").trim()
  return text ? text : fallback
}

function formatTaxonomy(taxonomy) {
  const raw = String(taxonomy ?? "").trim()
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
      String(value ?? "")
        .toLowerCase()
        .trim(),
    )
    .filter(Boolean)
    .join(" ")
}

function buildTooltipHTML(id, nodeMeta) {
  const organism = normalizeText(
    nodeMeta.organism || nodeMeta.species || id,
    id,
  )
  const genome = normalizeText(nodeMeta.genome)
  const group = normalizeText(nodeMeta.group)
  const environments = normalizeText(
    nodeMeta.environments || nodeMeta.environment,
  )
  const regulon = normalizeText(nodeMeta.regulon)
  const taxonomy = formatTaxonomy(nodeMeta.taxonomy)
  const operonHtml = SHOW_OPERON_IN_TOOLTIP ? renderOperonHTML(nodeMeta) : ""

  return `
    <div class="tooltip-card">
      <div class="tooltip-title">${escapeHtml(organism)}</div>
      <div class="tooltip-subtitle">${escapeHtml(id)}</div>

      <div class="tooltip-meta-grid">
        <div class="tooltip-key">Genome</div>
        <div class="tooltip-value">${escapeHtml(genome)}</div>

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

// load SVG + metadata + leaf labels + add tooltip interactivity
Promise.all([
  fetch("./images/tree.svg").then((r) => r.text()),
  fetch("./metadata.json").then((r) => r.json()),
])
  .then(([svg, metadata]) => {
    document.getElementById("tree-container").innerHTML = svg
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
    })

    console.log("metadata-linked labels:", labelData.length)
  })
  .catch(console.error)

// add search functionality

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
