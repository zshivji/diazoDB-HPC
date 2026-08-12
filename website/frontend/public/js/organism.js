const detailRoot = document.getElementById("organism-detail")

const GENE_COLORS = {
  nifA: "#c084fc",
  nifB: "#f472b6",
  nifH: "#60a5fa",
  nifD: "#f59e0b",
  nifK: "#10b981",
  nifE: "#f97316",
  nifN: "#14b8a6",
  nifV: "#a3e635",
  anfG: "#818cf8",
  anfO: "#fb7185",
  vnfG: "#22d3ee",
}

const DEFAULT_GENE_COLOR = "#9ca3af"

function escapeHtml(str) {
  return String(str ?? "")
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;")
    .replace(/'/g, "&#039;")
}

function normalizeText(value, fallback = "Not available") {
  const text = String(value ?? "").trim()
  return text ? text : fallback
}

function normalizeTreeIdentifier(value) {
  const normalized = String(value || "")
    .trim()
    .replace(/^'+|'+$/g, "")
    .replace(/\b(GB|RS)\s+(GCA|GCF)\s+(\d+\.\d+)\b/gi, "$1_$2_$3")
    .replace(/\b([A-Z]{2})\s+([A-Z]*\d+(?:\.\d+)?)\b/g, "$1_$2")

  return normalized
    .split("|")
    .map((part) => part.trim())
    .filter(Boolean)
    .join(" | ")
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

function getContigFromKey(id) {
  const parts = String(id || "").split("|")
  return normalizeText(parts[3] || parts[parts.length - 1])
}

function getEnvironment(nodeMeta) {
  const environment = String(nodeMeta.environments || nodeMeta.environment || "").trim()
  return !environment || environment.toLowerCase() === "none"
    ? "Not Available"
    : environment
}

function findOrganismRecords(metadata, requestedId) {
  const decodedId = decodeURIComponent(requestedId)
  const normalizedId = normalizeTreeIdentifier(decodedId)
  let selectedId = Object.prototype.hasOwnProperty.call(metadata, decodedId)
    ? decodedId
    : Object.prototype.hasOwnProperty.call(metadata, normalizedId)
      ? normalizedId
      : null

  const requestedLower = normalizedId.toLowerCase()
  const entry = selectedId
    ? [selectedId, metadata[selectedId]]
    : Object.entries(metadata).find(([key, meta]) => {
    const normalizedKey = normalizeTreeIdentifier(key).toLowerCase()
    const genome = String(meta.genome || "").toLowerCase()
    return (
      normalizedKey === requestedLower ||
      genome === requestedLower ||
      genome.endsWith(`_${requestedLower}`) ||
      normalizedKey.includes(`|${requestedLower}|`) ||
      normalizedKey.includes(`_${requestedLower}|`) ||
      normalizedKey.endsWith(`|${requestedLower}`)
    )
  })

  if (!entry) return null

  const selectedMeta = entry[1]
  const organism = String(
    selectedMeta.organism || selectedMeta.species || entry[0].split("|")[0],
  ).trim().toLowerCase()

  return Object.entries(metadata)
    .filter(([key, meta]) => {
      const recordOrganism = String(
        meta.organism || meta.species || key.split("|")[0],
      ).trim().toLowerCase()
      return recordOrganism === organism
    })
    .map(([id, meta]) => ({ id, meta }))
}

async function fetchJson(paths) {
  let lastError = null

  for (const path of paths) {
    try {
      const response = await fetch(path)
      if (!response.ok) {
        throw new Error(`${response.status} ${response.statusText}`)
      }
      const text = await response.text()
      const trimmedText = text.trim()
      if (trimmedText.startsWith("<")) {
        throw new Error(`${path} returned HTML instead of JSON`)
      }
      return JSON.parse(trimmedText.replace(/\bNaN\b/g, "null"))
    } catch (error) {
      lastError = new Error(`Unable to load ${path}: ${error.message}`)
    }
  }

  throw lastError
}

function renderScaleHTML(regionStart, regionEnd) {
  const fractions = [0, 0.2, 0.4, 0.6, 0.8, 1]

  return fractions
    .map((fraction, i) => {
      const value = Math.round(regionStart + fraction * (regionEnd - regionStart))
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
}

function renderOperonHTML(nodeMeta) {
  if (!nodeMeta?.operon?.genes || nodeMeta.operon.genes.length === 0) {
    return '<p class="organism-muted">No operon data available.</p>'
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
  const trackWidth = 404
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

  return `
    <div class="organism-operon">
      <div class="tooltip-operon-track">
        ${geneHtml}
      </div>
      <div class="tooltip-operon-scale">
        <div class="tooltip-operon-scale-line"></div>
        ${renderScaleHTML(regionStart, regionEnd)}
      </div>
    </div>
  `
}

function renderGeneTable(nodeMeta) {
  const genes = nodeMeta?.operon?.genes
  if (!genes || genes.length === 0) return ""

  const rows = [...genes]
    .sort((a, b) => Number(a.start) - Number(b.start))
    .map((gene) => `
      <tr>
        <td>${escapeHtml(gene.gene_name || "Not available")}</td>
        <td>${escapeHtml(gene.gene_id || "Not available")}</td>
        <td>${escapeHtml(normalizeDirection(gene.direction))}</td>
        <td>${escapeHtml(gene.start ?? "Not available")}</td>
        <td>${escapeHtml(gene.end ?? "Not available")}</td>
        <td>${escapeHtml(gene.product || "Not available")}</td>
      </tr>
    `)
    .join("")

  return `
    <section class="organism-section">
      <h2>Operon genes</h2>
      <div class="organism-table-wrap">
        <table class="organism-gene-table">
          <thead>
            <tr>
              <th>Gene</th>
              <th>Gene ID</th>
              <th>Direction</th>
              <th>Start</th>
              <th>End</th>
              <th>Product</th>
            </tr>
          </thead>
          <tbody>${rows}</tbody>
        </table>
      </div>
    </section>
  `
}

function renderOrganismPage(records) {
  const firstRecord = records[0]
  const organism = normalizeText(
    firstRecord.meta.organism || firstRecord.meta.species || firstRecord.id,
    firstRecord.id,
  )
  const genomes = [...new Set(records.map(({ meta }) => normalizeText(meta.genome)))]
  const contigs = [...new Set(records.map(({ id }) => getContigFromKey(id)))]
  const environments = [...new Set(records.map(({ meta }) => getEnvironment(meta)))]
  const groups = [...new Set(records.map(({ meta }) => normalizeText(meta.group)))]
  const regulons = [...new Set(records.map(({ meta }) => normalizeText(meta.regulon)))]
  const taxonomies = [...new Set(records.map(({ meta }) => formatTaxonomy(meta.taxonomy)))]

  document.title = `${organism} - DiazoDB`

  return `
    <div class="organism-back-row">
      <a href="/">Back to phylogeny</a>
    </div>

    <section class="organism-hero">
      <p class="organism-eyebrow">${escapeHtml(groups.join(", "))}</p>
      <h1>${escapeHtml(organism)}</h1>
      <p class="organism-id">${records.length} operon${records.length === 1 ? "" : "s"}</p>
    </section>

    <section class="organism-summary" aria-label="Organism metadata">
      <div>
        <span>Genome</span>
        <strong>${escapeHtml(genomes.join(", "))}</strong>
      </div>
      <div>
        <span>Contig</span>
        <strong>${escapeHtml(contigs.join(", "))}</strong>
      </div>
      <div>
        <span>Environment</span>
        <strong>${escapeHtml(environments.join(", "))}</strong>
      </div>
      <div>
        <span>Regulon</span>
        <strong>${escapeHtml(regulons.join(", "))}</strong>
      </div>
    </section>

    <section class="organism-section">
      <h2>Taxonomy</h2>
      <p class="organism-taxonomy">${escapeHtml(taxonomies.join("; "))}</p>
    </section>

    ${records.map(({ id, meta }, index) => `
      <section class="organism-section">
        <h2>Operon ${index + 1}</h2>
        <p class="organism-id">${escapeHtml(id)}</p>
        ${renderOperonHTML(meta)}
      </section>
      ${renderGeneTable(meta)}
    `).join("")}
  `
}

async function initOrganismPage() {
  const params = new URLSearchParams(window.location.search)
  const requestedId = params.get("id")

  if (!requestedId) {
    detailRoot.innerHTML = `
      <div class="organism-empty">
        <h1>No organism selected</h1>
        <p>Choose a node from the phylogeny tree to open an organism page.</p>
        <a href="/">Back to phylogeny</a>
      </div>
    `
    return
  }

  try {
    const metadata = await fetchJson([
      "./results/metadata.json",
      "../results/metadata.json",
      "./metadata.json",
    ])
    const records = findOrganismRecords(metadata, requestedId)

    if (!records || records.length === 0) {
      detailRoot.innerHTML = `
        <div class="organism-empty">
          <h1>Organism not found</h1>
          <p>No metadata record matched <code>${escapeHtml(requestedId)}</code>.</p>
          <a href="/">Back to phylogeny</a>
        </div>
      `
      return
    }

    detailRoot.innerHTML = renderOrganismPage(records)
  } catch (error) {
    console.error(error)
    detailRoot.innerHTML = `
      <div class="organism-empty">
        <h1>Unable to load organism</h1>
        <p>Metadata could not be loaded for this organism page.</p>
        <a href="/">Back to phylogeny</a>
      </div>
    `
  }
}

initOrganismPage()
