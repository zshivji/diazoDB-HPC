/* ============================================================
   js/article.js — TOC scroll-spy for article pages
   ============================================================ */

document.addEventListener("DOMContentLoaded", () => {
  const tocLinks = document.querySelectorAll(".toc-item")
  if (!tocLinks.length) return

  // Build a map of id → toc link
  const sectionMap = new Map()
  tocLinks.forEach((link) => {
    const href = link.getAttribute("href")
    if (href?.startsWith("#")) {
      const id = href.slice(1)
      const el = document.getElementById(id)
      if (el) sectionMap.set(el, link)
    }
  })

  const sections = Array.from(sectionMap.keys())

  function onScroll() {
    const scrollY = window.scrollY + 100 // offset for sticky header

    let active = sections[0]
    for (const section of sections) {
      if (section.offsetTop <= scrollY) {
        active = section
      }
    }

    tocLinks.forEach((l) => l.classList.remove("active"))
    const activeLink = sectionMap.get(active)
    if (activeLink) activeLink.classList.add("active")
  }

  window.addEventListener("scroll", onScroll, { passive: true })
  onScroll() // run once on load
})
