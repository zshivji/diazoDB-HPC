/* ============================================================
   js/main.js — Search, tab filtering, misc interactions
   ============================================================ */

/* ---- Tab filtering ---- */
function setFilter(filter) {
  const cards = document.querySelectorAll(".article-card")
  const noResults = document.getElementById("no-results")
  let visible = 0

  cards.forEach((card) => {
    const match = filter === "all" || card.dataset.category === filter
    card.style.display = match ? "" : "none"
    if (match) visible++
  })

  if (noResults) noResults.hidden = visible > 0
}

document.addEventListener("DOMContentLoaded", () => {
  // Tab clicks
  document.querySelectorAll(".tab").forEach((tab) => {
    tab.addEventListener("click", () => {
      document
        .querySelectorAll(".tab")
        .forEach((t) => t.classList.remove("active"))
      tab.classList.add("active")
      setFilter(tab.dataset.filter)
      // clear search
      const si = document.getElementById("search-input")
      if (si) si.value = ""
    })
  })

  // Search
  const searchInput = document.getElementById("search-input")
  if (searchInput) {
    searchInput.addEventListener("input", () => {
      const q = searchInput.value.trim().toLowerCase()
      const cards = document.querySelectorAll(".article-card")
      const noResults = document.getElementById("no-results")
      let visible = 0

      if (q === "") {
        // restore active tab filter
        const activeTab = document.querySelector(".tab.active")
        const filter = activeTab ? activeTab.dataset.filter : "all"
        setFilter(filter)
        return
      }

      // reset tabs when searching
      document
        .querySelectorAll(".tab")
        .forEach((t) => t.classList.remove("active"))

      cards.forEach((card) => {
        const match =
          card.dataset.title.includes(q) || card.dataset.excerpt.includes(q)
        card.style.display = match ? "" : "none"
        if (match) visible++
      })

      if (noResults) noResults.hidden = visible > 0
    })
  }

  // Keyboard shortcut ⌘K / Ctrl+K to focus search
  document.addEventListener("keydown", (e) => {
    if ((e.metaKey || e.ctrlKey) && e.key === "k") {
      e.preventDefault()
      const si = document.getElementById("search-input")
      if (si) si.focus()
    }
  })
})
