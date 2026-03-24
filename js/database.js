// define variables
const tbody = document.getElementById("tbody");
const tpl = document.getElementById("rowTemplate");
const viewport = document.getElementById("tableViewport");
const topSpacer = document.getElementById("topSpacer");
const bottomSpacer = document.getElementById("bottomSpacer");

const ROW_HEIGHT = 40; // must match CSS .data-row height
const BUFFER = 20;

let visibleRows = [];
let currentSort = { key: null, dir: "asc" };

const table = document.getElementById("diazoDB-table");
if (!tbody || !tpl || !viewport || !topSpacer || !bottomSpacer || !table) {
  throw new Error("Missing required table elements in database.html");
}

const headers = Array.from(table.querySelectorAll("thead th[data-column-name]"));

const columnKeyMap = {
  Organism: "Organism",
  NitrogenaseSet: "NitrogenaseSet",
  GroupNo: "GroupNo",
  Env: "Env",
  GeneCluster: "GeneCluster",
  Regulon: "Regulon",
  PredGrowthTemp: "PredGrowthTemp",
  GenomeAcc: "GenomeAcc",
  ContigAcc: "ContigAcc",
  GeneAcc: "GeneAcc",
  GTDBPhylo: "GTDBPhylo",
  NCBIPyhlo: "NCBIPyhlo"
};

const rowCellMap = {
  ".col-Organism": "Organism",
  ".col-NitrogenaseSet": "NitrogenaseSet",
  ".col-GroupNo": "GroupNo",
  ".col-Env": "Env",
  ".col-GeneCluster": "GeneCluster",
  ".col-Regulon": "Regulon",
  ".col-PredGrowthTemp": "PredGrowthTemp",
  ".col-ContigAcc": "ContigAcc",
  ".col-GeneAcc": "GeneAcc",
  ".col-GTDBPhylo": "GTDBPhylo",
  ".col-NCBIPyhlo": "NCBIPyhlo"
};

// parse CSV text into array of objects, using first line as keys
function parseCSVToObjects(csvText) {
  const lines = csvText.trim().split(/\r?\n/);
  const csvHeaders = lines[0].split(",").map(h => h.trim());

  return lines.slice(1)
    .filter(line => line.trim().length > 0)
    .map(line => {
      const cols = line.split(",");
      const obj = {};
      csvHeaders.forEach((header, i) => (obj[header] = (cols[i] ?? "").trim()));
      return obj;
    });
}

// load final nitrogenase data, display in table
function buildRowNode(row, idx) {
  const node = tpl.content.cloneNode(true);
  const tr = node.querySelector("tr");

  // checkbox + label wiring
  const input = tr.querySelector('input[type="checkbox"]');
  const label = tr.querySelector("label");
  const checkboxId = `_r_${idx}`;

  input.dataset.id = row.GeneAcc || row.GenomeAcc || String(idx);
  input.id = checkboxId;

  label.htmlFor = checkboxId;
  label.setAttribute("aria-label", input.dataset.id);

  // accession link
  const genomeAcc = String(row.GenomeAcc || "");
  const genomeParts = genomeAcc.split("_");
  const genomeAccRaw = genomeParts.length >= 3
    ? `${genomeParts[1]}_${genomeParts[2]}`
    : genomeAcc;
  const genomeA = tr.querySelector(".col-GenomeAcc a.BqBnJ");

  if (genomeA) {
    genomeA.textContent = genomeAccRaw;
    genomeA.href = genomeAccRaw
      ? `https://gtdb.ecogenomic.org/genome?gid=${encodeURIComponent(genomeAccRaw)}`
      : "#";
    genomeA.setAttribute("translate", "no");
  }

  // Fill row text cells using a selector-to-key map to keep field wiring centralized.
  Object.entries(rowCellMap).forEach(([selector, key]) => {
    const cell = tr.querySelector(selector);
    if (cell) {
      cell.textContent = row[key] || "";
    }
  });

  return node;
}

// virtualization
function renderVirtualRows() {
  const total = visibleRows.length;
  const scrollTop = viewport.scrollTop;
  const viewportHeight = viewport.clientHeight;

  const start = Math.max(0, Math.floor(scrollTop / ROW_HEIGHT) - BUFFER);
  const end = Math.min(
    total,
    Math.ceil((scrollTop + viewportHeight) / ROW_HEIGHT) + BUFFER
  );

  const topHeight = start * ROW_HEIGHT;
  const bottomHeight = (total - end) * ROW_HEIGHT;

  topSpacer.firstElementChild.style.height = `${topHeight}px`;
  bottomSpacer.firstElementChild.style.height = `${bottomHeight}px`;

  tbody.querySelectorAll(".data-row").forEach(row => row.remove());

  const fragment = document.createDocumentFragment();
  for (let i = start; i < end; i++) {
    fragment.appendChild(buildRowNode(visibleRows[i], i + 1));
  }

  bottomSpacer.before(fragment);
}

function renderTable(rows) {
  visibleRows = [...rows];
  viewport.scrollTop = 0;
  renderVirtualRows();
}

function isNumberLike(v) {
  const s = String(v ?? "").trim();
  return s !== "" && !Number.isNaN(Number(s));
}

function compareValues(a, b) {
  const A = String(a ?? "").trim();
  const B = String(b ?? "").trim();

  if (isNumberLike(A) && isNumberLike(B)) {
    return Number(A) - Number(B);
  }

  return A.localeCompare(B, undefined, { numeric: true, sensitivity: "base" });
}

function sortBy(key) {
  if (currentSort.key === key) {
    currentSort.dir = currentSort.dir === "asc" ? "desc" : "asc";
  } else {
    currentSort.key = key;
    currentSort.dir = "asc";
  }

  headers.forEach(th => {
    th.classList.remove("data-table__header-cell--ascend", "data-table__header-cell--descend");
    const thKey = columnKeyMap[th.dataset.columnName] ?? th.dataset.columnName;
    if (thKey === key) {
      th.classList.add(
        currentSort.dir === "asc"
          ? "data-table__header-cell--ascend"
          : "data-table__header-cell--descend"
      );
    }
  });

  const sorted = [...visibleRows].sort((r1, r2) => {
    const cmp = compareValues(r1[key], r2[key]);
    return currentSort.dir === "asc" ? cmp : -cmp;
  });

  visibleRows = sorted;
  renderVirtualRows();
}

headers.forEach(th => {
  th.addEventListener("click", () => {
    const key = columnKeyMap[th.dataset.columnName] ?? th.dataset.columnName;
    sortBy(key);
  });
});

let rafPending = false;
viewport.addEventListener("scroll", () => {
  if (!rafPending) {
    rafPending = true;
    requestAnimationFrame(() => {
      renderVirtualRows();
      rafPending = false;
    });
  }
});

// load CSV, parse, and add rows to table
fetch("DiazoDB-data.csv")
  .then(r => {
    if (!r.ok) throw new Error(`HTTP error! status: ${r.status}`);
    return r.text();
  })
  .then(csvText => {
    const rows = parseCSVToObjects(csvText);
    console.log(`Parsed ${rows.length} rows from CSV.`);
    renderTable(rows);
  })
  .catch(error => console.error("Error loading CSV:", error));