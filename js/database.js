console.log("database.js is running");

// define variables
const tbody = document.getElementById("tbody");
const tpl = document.getElementById("rowTemplate");

// load final nitrogenase data, display in table
function addRow(row, idx) {
  const node = tpl.content.cloneNode(true);
  const tr = node.querySelector("tr");

  // checkbox + label wiring
  const input = tr.querySelector('input[type="checkbox"]');
  const label = tr.querySelector("label");
  const checkboxId = `_r_${idx}`;
  console.log("Adding row with accession: " + checkboxId);

  input.dataset.id = row.GeneAcc || row.GenomeAcc || String(idx);
  input.id = checkboxId;

  label.htmlFor = checkboxId;
  label.setAttribute('aria-label', input.dataset.id);

  // accession link
  const genomeAccRaw = row.GenomeAcc.split("_")[1] + "_" + row.GenomeAcc.split("_")[2] || "";
  const genomeA = tr.querySelector(".col-GenomeAcc a.BqBnJ");

  if (genomeA) {
    genomeA.textContent = genomeAccRaw;
    genomeA.href = genomeAccRaw
      ? `https://gtdb.ecogenomic.org/genome?gid=${encodeURIComponent(genomeAccRaw)}`
      : "#";
    genomeA.setAttribute("translate", "no");
  }

  console.log("Adding row with link: " + genomeAccRaw);

  // plain text columns
  tr.querySelector(".col-Organism").textContent = row.Organism || "";
  tr.querySelector(".col-NitrogenaseSet").textContent = row.NitrogenaseSet || "";
  tr.querySelector(".col-GroupNo").textContent = row.GroupNo || "";
  tr.querySelector(".col-Env").textContent = row.Env || "";
  tr.querySelector(".col-GeneCluster").textContent = row.GeneCluster || ""; 
  tr.querySelector(".col-Regulon").textContent = row.Regulon || "";
  tr.querySelector(".col-PredGrowthTemp").textContent = row.PredGrowthTemp || "";
  tr.querySelector(".col-ContigAcc").textContent = row.ContigAcc || "";
  tr.querySelector(".col-GeneAcc").textContent = row.GeneAcc || "";
  tr.querySelector(".col-GTDBPhylo").textContent = row.GTDBPhylo || "";
  tr.querySelector(".col-NCBIPyhlo").textContent = row.NCBIPyhlo || "";

  tbody.appendChild(node);

}

// parse CSV text into array of objects, using first line as keys
function parseCSVToObjects(csvText) {
  const lines = csvText.trim().split(/\r?\n/);
  const headers = lines[0].split(",").map(h => h.trim());

  return lines.slice(1)
    .filter(line => line.trim().length >0)
    .map(line => {
    const cols = line.split(",");
    const obj = {};
    headers.forEach((header, i) => (obj[header] = (cols[i] ??"").trim()));
    return obj;
  });
}

// load CSV, parse, and add rows to table
fetch("DiazoDB-data.csv")
      .then(r => {
        if (!r.ok) throw new Error(`HTTP error! status: ${r.status}`);
        return r.text();
      })
      .then(csvText => {
        const rows = parseCSVToObjects(csvText);
        console.log(`Parsed ${rows.length} rows from CSV.`);
        rows.forEach((row, idx) => {
          console.log(`Adding row ${idx} with GeneAcc: ${row.GeneAcc}`);
          addRow(row, idx+1)});
        
      })
      .catch(error => console.error("Error loading CSV:", error));