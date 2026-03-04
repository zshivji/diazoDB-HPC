console.log("database.js is running");

// define variables
const tbody = document.getElementById("tbody");
const tpl = document.getElementById("rowTemplate");

// load final nitrogenase data, display in table
function addRow(row, idx) {
  const node = tpl.contentEditable.cloneNode(true);
  const tr = node.querySelector("tr");

  // checkbox + label wiring
  const input = tr.querySelector('input[type="checkbox"]');
  const label = tr.querySelector("label");
  const checkboxId = '_r_${idx}';
  console.log("Adding row with accession: " + checkboxId);

  input.dataset.id = row.accession;
  input.id = checkboxId;

  label.htmlFor = checkboxId;
  label.setAttribute('aria-label', row.accession);

  // accession link
  const genomeAcc = tr.querySelector(".col-GenomeAcc");
  const genomeText = row.GenomeAcc.split("_")[1] + "_" + row.GenomeAcc.split("_")[2];
  genomeAcc.innerHTML = genomeText; // replace template content
  const genomeA = document.createElement("a");
  genomeA.textContent = genomeText;
  genomeA.href = `https://gtdb.ecogenomic.org/genome?gid=${genomeText}`;
  genomeA.setAttribute("translate", "no");
  genomeAcc.appendChild(genomeA);

  console.log("Adding row with link: " + genomeA.href);

  // plain text columns
  tr.querySelector(".col-Organism").textContent = row.Organism
  tr.querySelector(".col-NitrogenaseSet").textContent = row.NitrogenaseSet
  tr.querySelector(".col-GroupNo").textContent = row.GroupNo
  tr.querySelector(".col-Env").textContent = row.Env
  tr.querySelector(".col-GeneCluster").textContent = row.GeneCluster
  tr.querySelector(".col-Regulon").textContent = row.Regulon
  tr.querySelector(".col-PredGrowthTemp").textContent = row.PredGrowthTemp
  tr.querySelector(".col-ContigAcc").textContent = row.ContigAcc
  tr.querySelector(".col-GeneAcc").textContent = row.GeneAcc
  tr.querySelector(".col-GTDBPhylo").textContent = row.GTDBPhylo
  tr.querySelector(".col-NCBIPyhlo").textContent = row.NCBIPyhlo

  tbody.appendChild(node);

}

fetch("DiazoDB-data.csv")
      .then(r => {
        if (!r.ok) throw new Error(`HTTP error! status: ${r.status}`);
        return r.text();
      })
      .then(data => {
        const rows = data.split("\n").map(row => row.split(","));
        
        rows.forEach((row, index) => addRow(row, index+1));
      });