/* ============================================================
   js/articles.js — Article registry & homepage renderer
   ============================================================
   HOW TO ADD A NEW ARTICLE:
   1. Add an entry to ARTICLES array below.
   2. Create the folder: articles/your-slug/
   3. Copy article-template.html into it as index.html
   4. Fill in your content.
   5. Optionally add a thumbnail: images/thumbs/your-slug.jpg
   ============================================================ */

const ARTICLES = [
  {
    slug:     "test-slug",
    title:    "Title",
    excerpt:  "One sentence summary",
    category: "protein-structure",
    categoryLabel: "Protein structure",
    readTime: 12,
    updated:  "March 2026",
    featured: false,
    thumb:    null,   // e.g. "images/thumbs/mofe-protein.jpg"
    organisms: ["A. vinelandii", "K. pneumoniae"],
  },
  {
    slug:     "test-slug",
    title:    "Title",
    excerpt:  "One sentence summary",
    category: "mechanism",
    categoryLabel: "Mechanism",
    readTime: 12,
    updated:  "March 2026",
    featured: false,
    thumb:    null,   // e.g. "images/thumbs/mofe-protein.jpg"
    organisms: ["A. vinelandii", "K. pneumoniae"],
  },
  {
    slug:     "test-slug",
    title:    "Title",
    excerpt:  "One sentence summary",
    category: "regulation",
    categoryLabel: "Regulation",
    readTime: 12,
    updated:  "March 2026",
    featured: false,
    thumb:    null,   // e.g. "images/thumbs/mofe-protein.jpg"
    organisms: ["A. vinelandii", "K. pneumoniae"],
  },
  {
    slug:     "test-slug",
    title:    "Title",
    excerpt:  "One sentence summary",
    category: "genetics",
    categoryLabel: "Genetics",
    readTime: 12,
    updated:  "March 2026",
    featured: false,
    thumb:    null,   // e.g. "images/thumbs/mofe-protein.jpg"
    organisms: ["A. vinelandii", "K. pneumoniae"],
  },
  {
    slug:     "test-slug",
    title:    "Title",
    excerpt:  "One sentence summary",
    category: "ecology",
    categoryLabel: "Ecology",
    readTime: 12,
    updated:  "March 2026",
    featured: false,
    thumb:    null,   // e.g. "images/thumbs/mofe-protein.jpg"
    organisms: ["A. vinelandii", "K. pneumoniae"],
  }
];

/* ---- Derived data ---- */
const CATEGORIES = [
  { id: "protein-structure", label: "Protein structure" },
  { id: "mechanism",         label: "Mechanism" },
  { id: "regulation",        label: "Regulation" },
  { id: "genetics",          label: "Genetics" },
  { id: "cofactors",         label: "Cofactors" },
  { id: "ecology",           label: "Ecology" },
];

const KEY_RESIDUES = [
  "α-HisXXX","β-CysXXX"
];

/* ============================================================
   Rendering helpers
   ============================================================ */

function articleURL(slug) {
  return `articles/${slug}/index.html`;
}

function renderThumb(article) {
  if (article.thumb) {
    return `<img src="${article.thumb}" alt="${article.title}" loading="lazy" />`;
  }
  return `<span class="thumb-placeholder">${article.categoryLabel}</span>`;
}

function renderCard(article) {
  return `
    <div class="article-card" data-category="${article.category}" data-title="${article.title.toLowerCase()}" data-excerpt="${article.excerpt.toLowerCase()}"
         onclick="window.location='${articleURL(article.slug)}'">
      <div>
        <div class="article-tag">${article.categoryLabel}</div>
        <div class="article-title">${article.title}</div>
        <div class="article-excerpt">${article.excerpt}</div>
        <div class="article-meta">
          <span>${article.readTime} min read</span>
          <span class="meta-dot"></span>
          <span>${article.categoryLabel}</span>
          <span class="meta-dot"></span>
          <span>Updated ${article.updated}</span>
        </div>
      </div>
      <div class="article-thumb">${renderThumb(article)}</div>
    </div>`;
}

function renderFeatured(article) {
  return `
    <div class="featured-card" onclick="window.location='${articleURL(article.slug)}'">
      <div class="f-tag">${article.categoryLabel}</div>
      <div class="f-title">${article.title}</div>
      <div class="f-meta">${article.readTime} min read</div>
    </div>`;
}

function countByCategory() {
  return ARTICLES.reduce((acc, a) => {
    acc[a.category] = (acc[a.category] || 0) + 1;
    return acc;
  }, {});
}

/* ============================================================
   Mount everything on DOMContentLoaded
   ============================================================ */
document.addEventListener("DOMContentLoaded", () => {

  // --- Article feed
  const feed = document.getElementById("article-feed");
  if (feed) {
    feed.innerHTML = ARTICLES.map(renderCard).join("");
  }

  // --- Featured sidebar
  const featured = document.getElementById("featured-articles");
  if (featured) {
    featured.innerHTML = ARTICLES
      .filter(a => a.featured)
      .map(renderFeatured).join("");
  }

  // --- Category list with counts
  const catList = document.getElementById("category-list");
  if (catList) {
    const counts = countByCategory();
    catList.innerHTML = CATEGORIES.map(c => `
      <li class="cat-item" data-filter="${c.id}">
        <span class="cat-name">${c.label}</span>
        <span class="cat-count">${counts[c.id] || 0}</span>
      </li>`).join("");

    catList.querySelectorAll(".cat-item").forEach(item => {
      item.addEventListener("click", () => {
        const f = item.dataset.filter;
        setFilter(f);
        // sync tab UI
        document.querySelectorAll(".tab").forEach(t => {
          t.classList.toggle("active", t.dataset.filter === f);
        });
      });
    });
  }

  // --- Residue cloud
  const cloud = document.getElementById("residue-cloud");
  if (cloud) {
    cloud.innerHTML = KEY_RESIDUES
      .map(r => `<span class="residue-pill">${r}</span>`)
      .join(" ");
  }
});
