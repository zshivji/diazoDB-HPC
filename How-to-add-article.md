# DiazoDb — Nitrogenase Wiki

A static wiki for all things nitrogenase, hosted on GitHub Pages.

---

## How to add a new article

**Step 1 — Register it in `js/articles.js`**

Open `js/articles.js` and add an entry to the `ARTICLES` array:

```js
{
  slug:          "your-slug",          // must match folder name in articles/
  title:         "Your Article Title",
  excerpt:       "One or two sentence summary shown on the homepage card.",
  category:      "mechanism",          // see CATEGORIES list in articles.js
  categoryLabel: "Mechanism",
  readTime:      15,                   // estimated minutes
  updated:       "April 2026",
  featured:      false,                // true = appears in sidebar featured box
  thumb:         null,                 // or "images/thumbs/your-slug.jpg"
  organisms:     ["A. vinelandii"],    // shown in article right sidebar
},
```

**Step 2 — Create the article folder**

```bash
cp -r articles/article-template articles/your-slug
```

**Step 3 — Edit `articles/your-slug/index.html`**

Things to update (all marked with comments in the template):
- `<title>` and meta description tags
- Breadcrumb category link and label
- TOC links (`href="#section-id"`) — one per `<h2>`, sub-item per `<h3>`
- `class="article-eyebrow"` text
- `<h1>` article title
- Lead paragraph
- Meta bar: date and read time
- Body content (h2, h3, p, figures, tables, callouts)
- References
- Right sidebar: organisms, related articles, external links

**Step 4 — Add images (optional)**

- Homepage thumbnail: `images/thumbs/your-slug.jpg` (88×66 px or similar 4:3)
- Article figures: `images/figures/descriptive-name.png`
- Reference the thumbnail in `articles.js`: `thumb: "images/thumbs/your-slug.jpg"`
- Reference figures in the article HTML: `<img src="../../images/figures/name.png" />`

---

## Component reference

### Callout box (key findings, definitions, caveats)
```html
<div class="callout">
  <div class="callout-label">Key finding</div>
  <p>Your text here.</p>
</div>
```

### Residue table
```html
<table class="residue-table">
  <thead><tr><th>Residue</th><th>Subunit</th><th>Role</th></tr></thead>
  <tbody>
    <tr>
      <td><span class="r-pill">α-His195</span></td>
      <td>NifD (α)</td>
      <td>Role description.</td>
    </tr>
  </tbody>
</table>
```

### Figure block (with image)
```html
<div class="figure-block">
  <img src="../../images/figures/your-figure.png" alt="Alt text" />
  <p class="figure-caption"><strong>Figure 1.</strong> Caption text.</p>
</div>
```

### Figure placeholder (while writing)
```html
<div class="figure-block">
  <div class="figure-placeholder">[ Figure placeholder ]</div>
  <p class="figure-caption"><strong>Figure 1.</strong> Caption text.</p>
</div>
```

### Organism tag
```html
<span class="organism-tag">A. vinelandii</span>
```

### Residue pill (inline in text or tables)
```html
<span class="r-pill">α-His195</span>
```