import type { CSSProperties } from "react"
import { useEffect, useRef, useState } from "react"

const API_BASE = `${import.meta.env.VITE_API_URL}/api/v1`

const theme = {
  bg: "#f7f6f2",
  bgCard: "#ffffff",
  bgInput: "#fbfeff",
  border: "rgba(30, 30, 20, 0.12)",
  borderHover: "rgba(30, 30, 20, 0.35)",
  green: "#776c4d",
  greenDim: "#5f5437",
  greenGlow: "rgba(119,108,77,0.12)",
  text: "#1a1a18",
  textMuted: "#888780",
  textDim: "#4a4a47",
  textWhite: "#fbfeff",
  red: "#b33a3a",
  redDim: "#8c1f1f",
  amber: "#9a6500",
}

const styles: Record<string, CSSProperties> = {
  page: {
    minHeight: "100vh",
    backgroundColor: theme.bg,
    backgroundImage: "url('/images/header.png')",
    backgroundPosition: "top center",
    backgroundRepeat: "no-repeat",
    backgroundSize: "cover",
    color: theme.text,
    display: "flex",
    flexDirection: "column",
    fontFamily: "Lato, system-ui, sans-serif",
  },
  header: {
    borderBottom: `1px solid ${theme.textWhite}`,
    color: "rgba(255, 255, 255, 0.633)",
    padding: "0.5em 1em",
    display: "flex",
    alignItems: "center",
    alignContent: "center",
    flexWrap: "wrap",
    gap: "0.5ch",
    minHeight: "3rem",
    boxSizing: "border-box",
  },
  headerNav: {
    display: "flex",
    alignItems: "center",
    flexWrap: "wrap",
    gap: "0.5ch",
    padding: "0 0.5ch",
  },
  headerLink: {
    borderBottom: "1px solid transparent",
    color: theme.textWhite,
    display: "inline-block",
    fontSize: "16px",
    fontWeight: 600,
    margin: "0.1rem 0.25rem",
    textDecoration: "none",
  },
  headerLinkActive: {
    borderBottom: "1px solid transparent",
    color: theme.textWhite,
    display: "inline-block",
    fontSize: "16px",
    fontWeight: 600,
    margin: "0.1rem 0.25rem",
    textDecoration: "none",
  },
  logo: {
    alignItems: "center",
    color: theme.textWhite,
    display: "flex",
    fontSize: "20px",
    fontWeight: 600,
    gap: "8px",
    marginBottom: 0,
    textDecoration: "none",
    whiteSpace: "nowrap",
  },
  logoBadge: {
    background: theme.green,
    borderRadius: "4px",
    color: "#fff",
    fontSize: "16px",
    fontWeight: 500,
    letterSpacing: "0.06em",
    padding: "2px 6px",
  },
  main: {
    background: theme.bg,
    boxSizing: "border-box",
    flex: 1,
    margin: 0,
    padding: "1.5rem clamp(1rem, 3vw, 2.5rem)",
    width: "100%",
  },
  h1: {
    fontSize: "clamp(2rem, 4vw, 3.25rem)",
    fontWeight: 700,
    color: theme.text,
    marginBottom: "8px",
    lineHeight: 1.05,
  },
  lead: {
    fontSize: "15px",
    color: theme.textDim,
    lineHeight: 1.7,
    marginBottom: "36px",
    maxWidth: "640px",
  },
  card: {
    background: theme.bgCard,
    border: `1px solid ${theme.border}`,
    borderRadius: "6px",
    padding: "28px",
    marginBottom: "20px",
  },
  cardTitle: {
    fontSize: "12px",
    fontWeight: 700,
    letterSpacing: "0.08em",
    color: theme.green,
    textTransform: "uppercase",
    marginBottom: "20px",
  },
  label: {
    display: "block",
    fontSize: "12px",
    color: theme.textDim,
    marginBottom: "8px",
    fontWeight: 700,
  },
  input: {
    width: "100%",
    background: theme.bgInput,
    border: `1px solid ${theme.border}`,
    borderRadius: "4px",
    color: theme.text,
    fontFamily: "ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace",
    fontSize: "13px",
    padding: "10px 14px",
    marginBottom: "20px",
    boxSizing: "border-box",
    outline: "none",
  },
  textarea: {
    width: "100%",
    background: theme.bgInput,
    border: `1px solid ${theme.border}`,
    borderRadius: "4px",
    color: theme.text,
    fontFamily: "ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace",
    fontSize: "12px",
    padding: "12px 14px",
    marginBottom: "8px",
    boxSizing: "border-box",
    resize: "vertical",
    minHeight: "160px",
    outline: "none",
    lineHeight: 1.6,
  },
  divider: {
    display: "flex",
    alignItems: "center",
    gap: "12px",
    margin: "16px 0",
    color: theme.textMuted,
    fontSize: "11px",
  },
  dividerLine: {
    flex: 1,
    height: "1px",
    background: theme.border,
  },
  checkbox: {
    display: "flex",
    alignItems: "flex-start",
    gap: "10px",
    marginBottom: "20px",
    cursor: "pointer",
  },
  checkboxLabel: {
    fontSize: "13px",
    color: theme.textDim,
    lineHeight: 1.5,
  },
  btn: {
    background: theme.green,
    color: "#fff",
    border: "none",
    borderRadius: "4px",
    padding: "12px 28px",
    fontSize: "13px",
    fontWeight: 700,
    fontFamily: "inherit",
    letterSpacing: "0.05em",
    cursor: "pointer",
  },
  btnSecondary: {
    background: "transparent",
    color: theme.textDim,
    border: `1px solid ${theme.border}`,
    borderRadius: "4px",
    padding: "10px 20px",
    fontSize: "12px",
    fontFamily: "inherit",
    cursor: "pointer",
  },
  cite: {
    background: "rgba(119,108,77,0.08)",
    border: `1px solid rgba(119,108,77,0.28)`,
    borderRadius: "6px",
    padding: "14px 18px",
    fontSize: "12px",
    color: theme.textMuted,
    marginBottom: "28px",
    lineHeight: 1.6,
  },
  row: {
    display: "grid",
    gridTemplateColumns: "repeat(auto-fit, minmax(min(100%, 320px), 1fr))",
    gap: "20px",
    alignItems: "start",
    maxWidth: "1100px",
  },
  sidebar: {
    background: theme.bgCard,
    border: `1px solid ${theme.border}`,
    borderRadius: "6px",
    padding: "22px",
  },
  sidebarTitle: {
    fontSize: "11px",
    fontWeight: 700,
    letterSpacing: "0.15em",
    color: theme.textMuted,
    textTransform: "uppercase",
    marginBottom: "14px",
  },
  sidebarText: {
    fontSize: "12px",
    color: theme.textDim,
    lineHeight: 1.7,
    marginBottom: "12px",
  },
  statusBox: {
    textAlign: "center",
    padding: "60px 20px",
  },
  spinner: {
    width: "40px",
    height: "40px",
    border: `2px solid ${theme.border}`,
    borderTop: `2px solid ${theme.green}`,
    borderRadius: "50%",
    margin: "0 auto 24px",
    animation: "spin 1s linear infinite",
  },
  table: {
    width: "100%",
    borderCollapse: "collapse" as const,
    fontSize: "12px",
  },
  th: {
    textAlign: "left" as const,
    padding: "10px 14px",
    borderBottom: `1px solid ${theme.border}`,
    color: theme.textMuted,
    fontWeight: 700,
    letterSpacing: "0.08em",
    fontSize: "11px",
    textTransform: "uppercase" as const,
  },
  td: {
    padding: "10px 14px",
    borderBottom: `1px solid ${theme.border}`,
    color: theme.textDim,
    fontFamily: "'IBM Plex Mono', monospace",
  },
}

function alertStyle(type: "info" | "warn" | "error"): CSSProperties {
  return {
    background:
      type === "info"
        ? "rgba(119,108,77,0.08)"
        : type === "warn"
          ? "rgba(154,101,0,0.08)"
          : "rgba(179,58,58,0.08)",
    border: `1px solid ${type === "info" ? theme.greenDim : type === "warn" ? theme.amber : theme.red}`,
    borderRadius: "6px",
    padding: "14px 18px",
    fontSize: "13px",
    color:
      type === "info"
        ? theme.textDim
        : type === "warn"
          ? theme.amber
          : theme.red,
    marginBottom: "24px",
    lineHeight: 1.6,
  }
}

// function badgeStyle(type: "success" | "error"): CSSProperties {
//   return {
//     display: "inline-block",
//     padding: "2px 8px",
//     borderRadius: "3px",
//     fontSize: "11px",
//     fontWeight: 700,
//     background:
//       type === "success" ? "rgba(119,108,77,0.16)" : "rgba(179,58,58,0.14)",
//     color: type === "success" ? theme.greenDim : theme.red,
//   }
// }

// --- Types ---
interface JobData {
  jobId: string
  email: string
  // useProdigal: boolean
}

interface UploadPageProps {
  onSubmit: (data: JobData) => void
}

interface WaitingPageProps {
  jobId: string
  onSuccess: (resultFilename: string) => void
  onFailure: () => void
}

interface ResultsPageProps {
  jobId: string
  resultFilename: string
  onReset: () => void
}

interface FailurePageProps {
  jobId: string
  onReset: () => void
}

interface JobStatusResponse {
  id: string
  filename: string
  file_size_bytes: number | null
  status:
    | "created"
    | "uploading"
    | "transferring"
    | "ready"
    | "processing"
    | "complete"
    | "failed"
  result_filename?: string | null
}

interface ResultsData {
  headers: string[]
  rows: string[][]
}

function resultDownloadUrl(jobId: string, resultFilename: string) {
  return `${API_BASE}/classify/${jobId}/results/${encodeURIComponent(resultFilename)}`
}

function parseCsv(text: string): ResultsData {
  const rows: string[][] = []
  let row: string[] = []
  let field = ""
  let quoted = false

  for (let i = 0; i < text.length; i += 1) {
    const char = text[i]
    const next = text[i + 1]

    if (quoted) {
      if (char === '"' && next === '"') {
        field += '"'
        i += 1
      } else if (char === '"') {
        quoted = false
      } else {
        field += char
      }
    } else if (char === '"') {
      quoted = true
    } else if (char === ",") {
      row.push(field)
      field = ""
    } else if (char === "\n") {
      row.push(field)
      rows.push(row)
      row = []
      field = ""
    } else if (char !== "\r") {
      field += char
    }
  }

  if (field || row.length) {
    row.push(field)
    rows.push(row)
  }

  const [headers = [], ...bodyRows] = rows.filter((items) =>
    items.some((item) => item.trim()),
  )
  return { headers, rows: bodyRows }
}

function jobIdFromUrl() {
  return new URLSearchParams(window.location.search).get("job")
}

function setJobUrl(jobId: string | null) {
  const url = new URL(window.location.href)
  if (jobId) {
    url.searchParams.set("job", jobId)
  } else {
    url.searchParams.delete("job")
  }
  window.history.replaceState({}, "", `${url.pathname}${url.search}${url.hash}`)
}

// --- Citing banner ---
function CitingBanner() {
  return (
    <div style={styles.cite}>
      <strong style={{ color: theme.green }}>Please cite:</strong> If you use
      DiazoDB for research, please cite our Github.{" "}
      {/* <span style={{ color: theme.textDim }}>
        Predicting and classifying nitrogenase sequences from metagenomes using
        DiazoDB.
      </span> */}
    </div>
  )
}

// --- Upload page ---
function UploadPage({ onSubmit }: UploadPageProps) {
  const [email, setEmail] = useState("")
  const [sequences, setSequences] = useState("")
  const [fileName, setFileName] = useState<string | null>(null)
  // const [useProdigal, setUseProdigal] = useState(false)
  const [dragging, setDragging] = useState(false)
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState<string | null>(null)
  const fileRef = useRef<HTMLInputElement>(null)

  const handleFile = (file: File | undefined) => {
    if (!file) return
    setFileName(file.name)
    const reader = new FileReader()
    reader.onload = (e) => setSequences(e.target?.result as string)
    reader.readAsText(file)
  }

  const handleDrop = (e: React.DragEvent) => {
    e.preventDefault()
    setDragging(false)
    handleFile(e.dataTransfer.files[0])
  }

  const handleSubmit = async () => {
    if (!sequences.trim())
      return setError("Please provide sequences or upload a FASTA file.")
    if (!email.trim())
      return setError(
        "Please provide an email address for results notification.",
      )
    setError(null)
    setLoading(true)
    try {
      const fileSizeBytes = new Blob([sequences]).size
      const res = await fetch(`${API_BASE}/classify/`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          user_email: email,
          filename: fileName ?? "sequences.fasta",
          file_size_bytes: fileSizeBytes,
          // use_prodigal: useProdigal,
          sequences,
        }),
      })
      if (!res.ok) {
        const contentType = res.headers.get("content-type") || ""
        let detail = ""
        if (contentType.includes("application/json")) {
          const body = await res.json()
          if (body?.detail) {
            detail = typeof body.detail === "string" ? body.detail : JSON.stringify(body.detail)
          } else {
            detail = JSON.stringify(body)
          }
        } else {
          detail = await res.text()
        }
        const suffix = detail ? ` - ${detail}` : ""
        throw new Error(`HTTP ${res.status} ${res.statusText}${suffix}`)
      }
      const data = await res.json()
      onSubmit({ jobId: data.id, email })
    } catch (err) {
      const message = err instanceof Error ? err.message : "Unknown error"
      setError(`Could not submit job: ${message}`)
    } finally {
      setLoading(false)
    }
  }

  return (
    <div style={styles.main}>
      <h1 style={styles.h1}>Classify Nitrogenase Sequences</h1>
      <p style={styles.lead}>
        DiazoDB provides accurate prediction and classification of nitrogenase
        sequences from metagenomes.
      </p>
      <CitingBanner />

      <div style={styles.row}>
        <div>
          <div style={styles.card}>
            <div style={styles.cardTitle}>Submit sequences</div>
            {error && <div style={alertStyle("error")}>{error}</div>}
            <div style={alertStyle("info")}>
              DiazoDB annotates nitrogenase catalytic proteins (NifH, NifD,
              NifK, NifE, NifK) and classifies them into Groups I-IV. 
            </div>

            <label style={styles.label}>
              Email address (for results notification)
            </label>
            <input
              style={styles.input}
              type="email"
              placeholder="you@institution.edu"
              value={email}
              onChange={(e) => setEmail(e.target.value)}
            />

            <label style={styles.label}>Paste FASTA sequences</label>
            <textarea
              style={styles.textarea}
              placeholder={
                ">seq1\nMRQCIAIGYKGGIGKSTTANTLAAMLAQHGAKVGLTGD...\n\n>seq2\n..."
              }
              value={sequences}
              onChange={(e) => {
                setSequences(e.target.value)
                setFileName(null)
              }}
            />

            <div style={styles.divider}>
              <div style={styles.dividerLine} />
              or upload a file
              <div style={styles.dividerLine} />
            </div>

            <div
              style={{
                border: `1px dashed ${dragging ? theme.green : theme.borderHover}`,
                borderRadius: "6px",
                padding: "24px",
                textAlign: "center",
                cursor: "pointer",
                marginBottom: "20px",
                background: dragging ? theme.greenGlow : theme.bgInput,
              }}
              onClick={() => fileRef.current?.click()}
              onDragOver={(e) => {
                e.preventDefault()
                setDragging(true)
              }}
              onDragLeave={() => setDragging(false)}
              onDrop={handleDrop}
            >
              <input
                ref={fileRef}
                type="file"
                accept=".fasta,.fa,.fna,.faa,.txt"
                style={{ display: "none" }}
                onChange={(e) => handleFile(e.target.files?.[0])}
              />
              {fileName ? (
                <span style={{ color: theme.green, fontSize: "13px" }}>
                  ✓ {fileName}
                </span>
              ) : (
                <span style={{ color: theme.textMuted, fontSize: "12px" }}>
                  Drop a FASTA file here, or click to browse
                </span>
              )}
            </div>

            {/* <div style={styles.divider}>
              <div style={styles.dividerLine} />
              options
              <div style={styles.dividerLine} />
            </div> */}

            {/* <label style={styles.checkbox}>
              <input
                type="checkbox"
                checked={useProdigal}
                onChange={() => setUseProdigal(!useProdigal)}
                style={{ marginTop: "2px", accentColor: theme.green }}
              />
              <div>
                <div
                  style={{
                    ...styles.checkboxLabel,
                    color: theme.text,
                    fontWeight: 600,
                  }}
                >
                  Predict CDS with Prodigal
                </div>
                <div style={styles.checkboxLabel}>
                  If your input is nucleotide (DNA) sequences, DiazoDB can run
                  Prodigal to predict coding sequences before classification. Do
                  not check this for protein sequences.
                </div>
              </div>
            </label> */}

            <button
              style={{ ...styles.btn, opacity: loading ? 0.6 : 1 }}
              onClick={handleSubmit}
              disabled={loading}
            >
              {loading ? "Submitting..." : "Submit →"}
            </button>
          </div>
        </div>

        <div>
          <div style={styles.sidebar}>
            <div style={styles.sidebarTitle}>Instructions</div>
            <p style={styles.sidebarText}>
              Paste FASTA-formatted sequences or upload a <code>.fasta</code>{" "}
              file.
            </p>
            {/* <p style={styles.sidebarText}>
              For nucleotide input from metagenomes, enable{" "}
              <strong style={{ color: theme.textDim }}>
                Prodigal CDS prediction
              </strong>{" "}
              to automatically predict protein-coding sequences before
              classification.
            </p> */}
            <p style={styles.sidebarText}>
              Provide your email to receive a link when the job completes or
              fails.
            </p>
          </div>

          <div style={{ ...styles.sidebar, marginTop: "16px" }}>
            <div style={styles.sidebarTitle}>Limits</div>
            <p style={styles.sidebarText}>Max runtime: 72 hours</p>
            <p style={styles.sidebarText}>Max file size: 20GB</p>
            <p style={styles.sidebarText}>Results stored for: 2 weeks</p>
            <p
              style={{
                ...styles.sidebarText,
                color: theme.textMuted,
                fontSize: "11px",
              }}
            >
              Download results promptly — they will be delted after two weeks or in the event of a
              server crash.
            </p>
          </div>
        </div>
      </div>
    </div>
  )
}

// --- Waiting page ---
function WaitingPage({ jobId, onSuccess, onFailure }: WaitingPageProps) {
  const [lastUpdated, setLastUpdated] = useState<string | null>(null)
  const [elapsed, setElapsed] = useState(0)
  const [jobStatus, setJobStatus] = useState("checking")

  useEffect(() => {
    const start = Date.now()
    const timer = setInterval(
      () => setElapsed(Math.floor((Date.now() - start) / 1000)),
      1000,
    )
    const checkStatus = async () => {
      try {
        const res = await fetch(`${API_BASE}/classify/${jobId}`)
        if (res.ok) {
          setLastUpdated(new Date().toLocaleTimeString())
          const data = (await res.json()) as JobStatusResponse
          setJobStatus(data.status)
          if (data.status === "complete") {
            onSuccess(data.result_filename ?? "nif_clusters.csv")
          }
          if (data.status === "failed") onFailure()
        } else if (res.status === 404) {
          setJobStatus("not found")
        }
      } catch {}
    }
    void checkStatus()
    const poll = setInterval(checkStatus, 3000)
    return () => {
      clearInterval(timer)
      clearInterval(poll)
    }
  }, [jobId])

  const mins = Math.floor(elapsed / 60)
  const secs = elapsed % 60

  return (
    <div style={styles.main}>
      <style>{`@keyframes spin { to { transform: rotate(360deg); } }`}</style>
      <CitingBanner />
      <div style={{ ...styles.card, ...styles.statusBox }}>
        <div style={styles.spinner} />
        <div style={styles.cardTitle}>Processing</div>
        <div
          style={{
            fontSize: "11px",
            color: theme.textMuted,
            letterSpacing: "0.1em",
            marginBottom: "8px",
          }}
        >
          Job ID
        </div>
        <div
          style={{
            fontSize: "13px",
            color: theme.green,
            fontFamily: "monospace",
            marginBottom: "24px",
            letterSpacing: "0.05em",
          }}
        >
          {jobId}
        </div>

        <div
          style={{
            display: "grid",
            gridTemplateColumns: "1fr 1fr 1fr",
            gap: "16px",
            maxWidth: "480px",
            margin: "0 auto 32px",
          }}
        >
          {(
            [
              ["Elapsed", `${mins}m ${secs < 10 ? "0" : ""}${secs}s`],
              ["Last checked", lastUpdated ?? "—"],
              ["Status", jobStatus],
            ] as [string, string][]
          ).map(([k, v]) => (
            <div
              key={k}
              style={{
                background: theme.bgInput,
                border: `1px solid ${theme.border}`,
                borderRadius: "6px",
                padding: "14px 10px",
              }}
            >
              <div
                style={{
                  fontSize: "10px",
                  color: theme.textMuted,
                  letterSpacing: "0.1em",
                  textTransform: "uppercase",
                  marginBottom: "6px",
                }}
              >
                {k}
              </div>
              <div style={{ fontSize: "13px", color: theme.textDim }}>{v}</div>
            </div>
          ))}
        </div>

        <p
          style={{
            fontSize: "12px",
            color: theme.textMuted,
            marginBottom: "20px",
          }}
        >
          This page refreshes automatically. You can bookmark this URL and
          return any time.
        </p>
      </div>

    </div>
  )
}

// --- Results page ---
function ResultsPage({ jobId, resultFilename, onReset }: ResultsPageProps) {
  const [results, setResults] = useState<ResultsData | null>(null)
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)

  useEffect(() => {
    const loadResults = async () => {
      setLoading(true)
      setError(null)
      try {
        const res = await fetch(resultDownloadUrl(jobId, resultFilename))
        if (!res.ok) {
          throw new Error(`HTTP ${res.status} ${res.statusText}`)
        }
        const text = await res.text()
        setResults(parseCsv(text))
      } catch (err) {
        const message = err instanceof Error ? err.message : "Unknown error"
        setError(`Could not load results: ${message}`)
      } finally {
        setLoading(false)
      }
    }

    void loadResults()
  }, [jobId, resultFilename])

  const visibleHeaders = results?.headers ?? []

  return (
    <div style={styles.main}>
      <CitingBanner />
      <div
        style={{
          display: "flex",
          alignItems: "baseline",
          gap: "14px",
          marginBottom: "28px",
        }}
      >
        <h1 style={{ ...styles.h1, marginBottom: 0 }}>Results</h1>
        <span
          style={{
            fontSize: "12px",
            color: theme.textMuted,
            fontFamily: "monospace",
          }}
        >
          {jobId}
        </span>
      </div>

      <div style={styles.card}>
        <div
          style={{
            display: "flex",
            justifyContent: "space-between",
            alignItems: "center",
            marginBottom: "20px",
          }}
        >
          <div style={styles.cardTitle}>Classification Results</div>
          <button
            style={styles.btn}
            onClick={() =>
              window.open(resultDownloadUrl(jobId, resultFilename), "_blank")}
          >
            Download CSV ↓
          </button>
        </div>

        {loading && <div style={alertStyle("info")}>Loading results...</div>}
        {error && <div style={alertStyle("error")}>{error}</div>}
        {!loading && !error && results && results.rows.length === 0 && (
          <div style={alertStyle("warn")}>
            The job completed, but the result file did not contain any rows.
          </div>
        )}
        {!loading && !error && results && results.rows.length > 0 && (
          <>
            <div style={{ ...styles.cardTitle, color: theme.textMuted }}>
              Classification output
            </div>
            <div style={{ overflowX: "auto" }}>
              <table style={styles.table}>
                <thead>
                  <tr>
                    {visibleHeaders.map((header, index) => (
                      <th key={`${header}-${index}`} style={styles.th}>
                        {header || `Column ${index + 1}`}
                      </th>
                    ))}
                  </tr>
                </thead>
                <tbody>
                  {results.rows.map((row, rowIndex) => (
                    <tr key={`${row[0] ?? "row"}-${rowIndex}`}>
                      {visibleHeaders.map((_, cellIndex) => (
                        <td key={cellIndex} style={styles.td}>
                          {row[cellIndex] ?? ""}
                        </td>
                      ))}
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </>
        )}
      </div>

      <button style={styles.btnSecondary} onClick={onReset}>
        ← Submit another job
      </button>
    </div>
  )
}

// --- Failure page ---
function FailurePage({ jobId, onReset }: FailurePageProps) {
  return (
    <div style={styles.main}>
      <div
        style={{
          display: "flex",
          alignItems: "baseline",
          gap: "14px",
          marginBottom: "28px",
        }}
      >
        <h1 style={{ ...styles.h1, marginBottom: 0, color: theme.red }}>
          Job failed
        </h1>
        <span
          style={{
            fontSize: "12px",
            color: theme.textMuted,
            fontFamily: "monospace",
          }}
        >
          {jobId}
        </span>
      </div>

      <div style={alertStyle("error")}>
        An error occurred while running your job. Please check that your input
        is valid FASTA format. If the error persists, contact{" "}
        <a
          href={`mailto:admin@diazo.db?subject=DiazoDB Error: ${jobId}`}
          style={{ color: theme.red }}
        >
          the administrator
        </a>{" "}
        with your job ID.
      </div>

      <div style={styles.card}>
        <div style={styles.cardTitle}>Common causes</div>
        {[
          [
            "Invalid FASTA format",
            "Each sequence must start with > followed by a sequence ID on a new line.",
          ],
          /* [
            "Non-protein sequences without Prodigal",
            "If submitting DNA/nucleotide sequences, enable the Prodigal CDS prediction option.",
          ], */
          [
            "Empty or corrupt file",
            "Try re-exporting your FASTA file and resubmitting.",
          ],
        ].map(([title, desc]) => (
          <div
            key={title}
            style={{
              marginBottom: "18px",
              paddingBottom: "18px",
              borderBottom: `1px solid ${theme.border}`,
            }}
          >
            <div
              style={{
                fontSize: "13px",
                color: theme.text,
                marginBottom: "4px",
              }}
            >
              {title}
            </div>
            <div style={{ fontSize: "12px", color: theme.textMuted }}>
              {desc}
            </div>
          </div>
        ))}
      </div>

      <button style={styles.btn} onClick={onReset}>
        ← Try again
      </button>
    </div>
  )
}

// --- App shell ---
type View = "upload" | "waiting" | "results" | "failure"


export default function DiazoDB() {
  useEffect(() => {
    const measurementId = "G-89SKEBYC49"
    const scriptId = "google-gtag"
    const dataLayerWindow = window as Window & {
      dataLayer?: unknown[][]
    }

    dataLayerWindow.dataLayer = dataLayerWindow.dataLayer || []
    const gtag = (...args: unknown[]) => {
      dataLayerWindow.dataLayer?.push(args)
    }
    gtag("js", new Date())
    gtag("config", measurementId)

    if (!document.getElementById(scriptId)) {
      const script = document.createElement("script")
      script.id = scriptId
      script.async = true
      script.src = `https://www.googletagmanager.com/gtag/js?id=${measurementId}`
      document.head.appendChild(script)
    }
  }, [])

  const [initialJobId] = useState(() => jobIdFromUrl())
  const [view, setView] = useState<View>(() =>
    initialJobId ? "waiting" : "upload",
  )
  const [jobData, setJobData] = useState<JobData | null>(() =>
    initialJobId ? { jobId: initialJobId, email: "" } : null,
  )
  const [resultFilename, setResultFilename] = useState("nif_clusters.csv")

  const handleSubmit = (data: JobData) => {
    setJobData(data)
    setResultFilename("nif_clusters.csv")
    setJobUrl(data.jobId)
    setView("waiting")
  }

  const handleReset = () => {
    setJobData(null)
    setResultFilename("nif_clusters.csv")
    setJobUrl(null)
    setView("upload")
  }

  return (
    <div style={styles.page}>
      <header style={styles.header}>
        <a href="/" style={styles.logo}>
          DiazoDB <span style={styles.logoBadge}>Upload</span>
        </a>
        <nav style={styles.headerNav} aria-label="Primary navigation">
          <a href="/" style={styles.headerLink}>
            Phylogeny
          </a>
          <a href="/database" style={styles.headerLink}>
            DiazoDB
          </a>
          <a href="/classify" style={styles.headerLinkActive}>
            Upload
          </a>
          {/* <a href="/wiki" style={styles.headerLink}>Nitrogenase Wiki</a> */}
        </nav>
      </header>

      {view === "upload" && <UploadPage onSubmit={handleSubmit} />}
      {view === "waiting" && (
        <WaitingPage
          jobId={jobData?.jobId ?? ""}
          onSuccess={(filename) => {
            setResultFilename(filename)
            setView("results")
          }}
          onFailure={() => setView("failure")}
        />
      )}
      {view === "results" && (
        <ResultsPage
          jobId={jobData?.jobId ?? ""}
          resultFilename={resultFilename}
          onReset={handleReset}
        />
      )}
      {view === "failure" && (
        <FailurePage
          jobId={jobData?.jobId ?? ""}
          onReset={handleReset}
        />
      )}
    </div>
  )
}
