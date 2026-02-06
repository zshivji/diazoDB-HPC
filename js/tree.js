console.log("tree.js is running");

// grab tooltip element
const tooltip = document.getElementById("tooltip");

// load SVG + metadata + leaf labels + add tooltip interactivity
Promise.all([
  fetch("./tree.svg").then(r => r.text()),
  fetch("./metadata.json").then(r => r.json())
]).then(([svg, metadata]) => {
  document.getElementById("tree-container").innerHTML = svg;
  console.log("SVG + metadata loaded");

  // grab iTOL leaf labels (<text> elements)
  const labels = document.querySelectorAll("svg text");
  console.log("labels:", labels.length);

  const labelData = [];

  labels.forEach(label => {
    const id = label.textContent.trim().replace(/'/g, ""); // remove quotes
    console.log("SVG label:", id);
    console.log("metadata has key?", metadata[id]);
    
    if (!metadata[id]) return; // skip label if not in metadata

    labelData.push({
      label,
      id,
      text: [
        id,
        metadata[id].species,
        metadata[id].taxonomy,
        metadata[id].environment
      ].join("<br>").toLowerCase()
    });

    label.classList.add("tree-label");

    label.addEventListener("mouseover", e => {
      tooltip.innerHTML = `
        <b>${id}</b><br>
        ${metadata[id].species}<br>
        ${metadata[id].taxonomy}<br>
        ${metadata[id].environment}
      `;
      tooltip.style.display = "block";
    });

    label.addEventListener("mousemove", e => {
      tooltip.style.left = e.pageX + 10 + "px";
      tooltip.style.top = e.pageY + 10 + "px";
    });

    label.addEventListener("mouseout", () => {
      tooltip.style.display = "none";
    });
  });
}).catch(console.error);

// add search functionality
const searchInput = document.createElement("tree-search");

searchInput.addEventListener("input", () => {
  const query = searchInput.value.toLowerCase().trim();

  if (!query) {
    labelData.forEach(d => {
      d.label.classList.remove("tree-match", "tree-dim");
    });
    return;
  }

  labelData.forEach(d => {
    if (d.text.includes(query)) {
      d.label.classList.add("tree-match");
      d.label.classList.remove("tree-dim");
    } else {
      d.label.classList.remove("tree-match");
      d.label.classList.add("tree-dim");
    }
  });
});