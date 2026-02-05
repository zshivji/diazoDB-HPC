console.log("tree.js is running");


// // grab iTOL leaf labels (<text> elements)
// function attachTooltips(metadata) {
//   const labels = document.querySelectorAll("svg text");
//   console.log("labels:", labels.length);
// }

const tooltip = document.getElementById("tooltip");

Promise.all([
  fetch("./tree.svg").then(r => r.text()),
  fetch("./metadata.json").then(r => r.json())
]).then(([svg, metadata]) => {
  document.getElementById("tree-container").innerHTML = svg;
  console.log("SVG + metadata loaded");

  // grab iTOL leaf labels (<text> elements)
  const labels = document.querySelectorAll("svg text");
  console.log("labels:", labels.length);

  labels.forEach(label => {
    const id = label.textContent.trim(); // remove whitespace
    console.log("SVG label:", id);
    console.log("metadata has key?", metadata[id]);
    
    if (!metadata[id]) return; // skip label if not in metadata

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


