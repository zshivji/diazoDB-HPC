console.log("tree.js is running");

fetch("./tree.svg")
  .then(r => r.text())
  .then(svg => {
    document.getElementById("tree-container").innerHTML = svg;
    console.log("SVG injected");
  })
  .catch(console.error);

// const tooltip = document.getElementById("tooltip");

// Promise.all([
//   fetch("./tree.svg").then(r => r.text()).then(svgText => {
//     document.getElementById("tree-container").innerHTML = svgText;
//   }),
//   fetch("./metadata.json").then(r => r.json())
// ]).then(([svgText, metadata]) => {

//   document.getElementById("tree-container").innerHTML = svgText;

//   // iTOL leaf labels are <text> elements
//   const labels = document.querySelectorAll("svg text");

//   labels.forEach(label => {
//     const id = label.textContent.trim();

//     if (!metadata[id]) return;

//     label.classList.add("tree-label");

//     label.addEventListener("mouseover", e => {
//       tooltip.innerHTML = `
//         <b>${id}</b><br>
//         ${metadata[id].species}<br>
//         ${metadata[id].taxonomy}<br>
//         ${metadata[id].environment}
//       `;
//       tooltip.style.display = "block";
//     });

//     label.addEventListener("mousemove", e => {
//       tooltip.style.left = e.pageX + 10 + "px";
//       tooltip.style.top = e.pageY + 10 + "px";
//     });

//     label.addEventListener("mouseout", () => {
//       tooltip.style.display = "none";
//     });
//   });
// });
