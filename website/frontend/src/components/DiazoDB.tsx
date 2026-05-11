import type { CSSProperties } from "react"
import { useEffect, useRef, useState } from "react"

const API_BASE = `${import.meta.env.VITE_API_URL}/api/v1`

const theme = {
  bg: "#0a0e0f",
  bgCard: "#111618",
  bgInput: "#0d1214",
  border: "#1e2c2f",
  borderHover: "#2e4a4f",
  green: "#3dffa0",
  greenDim: "#1a8c56",
  greenGlow: "rgba(61,255,160,0.15)",
  text: "#e8f0ed",
  textMuted: "#5a7a6e",
  textDim: "#8aa89e",
  red: "#ff5c5c",
  redDim: "#8c1f1f",
  amber: "#ffb347",
}

const styles: Record<string, CSSProperties> = {
  page: {
    minHeight: "100vh",
    background: theme.bg,
    color: theme.text,
    fontFamily: "'IBM Plex Mono', 'Courier New', monospace",
  },
  header: {
    borderBottom: `1px solid ${theme.border}`,
    padding: "20px 40px",
    display: "flex",
    alignItems: "baseline",
    gap: "16px",
  },
  logo: {
    fontSize: "22px",
    fontWeight: 700,
    color: theme.green,
    letterSpacing: "-0.5px",
  },
  logoSub: {
    fontSize: "12px",
    color: theme.textMuted,
    letterSpacing: "0.08em",
  },
  main: {
    maxWidth: "860px",
    margin: "0 auto",
    padding: "48px 24px",
  },
  h1: {
    fontSize: "28px",
    fontWeight: 700,
    color: theme.text,
    marginBottom: "8px",
    letterSpacing: "-0.5px",
  },
  lead: {
    fontSize: "14px",
    color: theme.textDim,
    lineHeight: 1.7,
    marginBottom: "36px",
    maxWidth: "640px",
  },
  card: {
    background: theme.bgCard,
    border: `1px solid ${theme.border}`,
    borderRadius: "8px",
    padding: "28px",
    marginBottom: "20px",
  },
  cardTitle: {
    fontSize: "11px",
    fontWeight: 700,
    letterSpacing: "0.15em",
    color: theme.green,
    textTransform: "uppercase",
    marginBottom: "20px",
  },
  label: {
    display: "block",
    fontSize: "12px",
    color: theme.textDim,
    marginBottom: "8px",
    letterSpacing: "0.05em",
  },
  input: {
    width: "100%",
    background: theme.bgInput,
    border: `1px solid ${theme.border}`,
    borderRadius: "4px",
    color: theme.text,
    fontFamily: "inherit",
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
    fontFamily: "inherit",
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
    color: "#0a0e0f",
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
    background: "rgba(61,255,160,0.04)",
    border: `1px solid ${theme.greenDim}`,
    borderRadius: "6px",
    padding: "14px 18px",
    fontSize: "12px",
    color: theme.textMuted,
    marginBottom: "28px",
    lineHeight: 1.6,
  },
  row: {
    display: "grid",
    gridTemplateColumns: "1fr 340px",
    gap: "20px",
    alignItems: "start",
  },
  sidebar: {
    background: theme.bgCard,
    border: `1px solid ${theme.border}`,
    borderRadius: "8px",
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
        ? "rgba(61,255,160,0.06)"
        : type === "warn"
          ? "rgba(255,179,71,0.08)"
          : "rgba(255,92,92,0.08)",
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

function badgeStyle(type: "success" | "error"): CSSProperties {
  return {
    display: "inline-block",
    padding: "2px 8px",
    borderRadius: "3px",
    fontSize: "11px",
    fontWeight: 700,
    background:
      type === "success" ? "rgba(61,255,160,0.15)" : "rgba(255,92,92,0.15)",
    color: type === "success" ? theme.green : theme.red,
  }
}

// --- Types ---
interface JobData {
  jobId: string
  email: string
  useProdigal: boolean
}

interface UploadPageProps {
  onSubmit: (data: JobData) => void
}

interface WaitingPageProps {
  jobId: string
  onSuccess: () => void
  onFailure: () => void
}

interface ResultsPageProps {
  jobId: string
  onReset: () => void
}

interface FailurePageProps {
  jobId: string
  onReset: () => void
}

// --- Citing banner ---
function CitingBanner() {
  return (
    <div style={styles.cite}>
      <strong style={{ color: theme.green }}>Please cite:</strong> If you use
      DiazoDB for research, please cite our paper.{" "}
      <span style={{ color: theme.textDim }}>
        Predicting and classifying nitrogenase sequences from metagenomes using
        DiazoDB.
      </span>
    </div>
  )
}

// --- Upload page ---
function UploadPage({ onSubmit }: UploadPageProps) {
  const [email, setEmail] = useState("")
  const [sequences, setSequences] = useState("")
  const [fileName, setFileName] = useState<string | null>(null)
  const [useProdigal, setUseProdigal] = useState(false)
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
          use_prodigal: useProdigal,
        }),
      })
      if (!res.ok) throw new Error("Submission failed")
      const data = await res.json()
      onSubmit({ jobId: data.id, email, useProdigal })
    } catch {
      setError("Could not submit job. Make sure the backend is running.")
    } finally {
      setLoading(false)
    }
  }

  return (
    <div style={styles.main}>
      <h1 style={styles.h1}>Classify Nitrogenase Sequences</h1>
      <p style={styles.lead}>
        DiazoDB provides accurate prediction and classification of nitrogenase
        sequences from metagenomes. Supports direct sequence input and CDS
        prediction from raw nucleotide data via Prodigal.
      </p>
      <CitingBanner />

      <div style={styles.row}>
        <div>
          <div style={styles.card}>
            <div style={styles.cardTitle}>Submit sequences</div>
            {error && <div style={alertStyle("error")}>{error}</div>}
            <div style={alertStyle("info")}>
              DiazoDB classifies catalytic subunits of nitrogenases (NifH, AnfH,
              VnfH). Sequences that do not encode nitrogenase catalytic subunits
              will still be classified but results may be unreliable.
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

            <div style={styles.divider}>
              <div style={styles.dividerLine} />
              options
              <div style={styles.dividerLine} />
            </div>

            <label style={styles.checkbox}>
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
            </label>

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
            <p style={styles.sidebarText}>
              For nucleotide input from metagenomes, enable{" "}
              <strong style={{ color: theme.textDim }}>
                Prodigal CDS prediction
              </strong>{" "}
              to automatically predict protein-coding sequences before
              classification.
            </p>
            <p style={styles.sidebarText}>
              Provide your email to receive a link when the job completes or
              fails.
            </p>
          </div>

          <div style={{ ...styles.sidebar, marginTop: "16px" }}>
            <div style={styles.sidebarTitle}>Limits</div>
            <p style={styles.sidebarText}>Max runtime: 2 hours</p>
            <p style={styles.sidebarText}>Max sequences: 4,000</p>
            <p style={styles.sidebarText}>Results stored for: 2 weeks</p>
            <p
              style={{
                ...styles.sidebarText,
                color: theme.textMuted,
                fontSize: "11px",
              }}
            >
              Download results promptly — they may be deleted in the event of a
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

  useEffect(() => {
    const start = Date.now()
    const timer = setInterval(
      () => setElapsed(Math.floor((Date.now() - start) / 1000)),
      1000,
    )
    const poll = setInterval(async () => {
      try {
        const res = await fetch(`${API_BASE}/classify/${jobId}`)
        if (res.ok) {
          setLastUpdated(new Date().toLocaleTimeString())
          // When you add a status field, check it here:
          // const data = await res.json()
          // if (data.status === "SUCCESS") onSuccess()
          // if (data.status === "FAILURE") onFailure()
        }
      } catch {}
    }, 3000)
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
              ["Status", "pending"],
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

      {/* Demo buttons — remove when you have a real status endpoint */}
      <div
        style={{
          textAlign: "center",
          marginTop: "16px",
          display: "flex",
          gap: "12px",
          justifyContent: "center",
        }}
      >
        <button style={styles.btnSecondary} onClick={onSuccess}>
          Simulate success →
        </button>
        <button
          style={{
            ...styles.btnSecondary,
            borderColor: theme.redDim,
            color: theme.red,
          }}
          onClick={onFailure}
        >
          Simulate failure →
        </button>
      </div>
    </div>
  )
}

// --- Results page ---
function ResultsPage({ jobId, onReset }: ResultsPageProps) {
  const mockResults = [
    {
      id: "seq1",
      prediction: "NifH Group I",
      confidence: "0.97",
      evalue: "1e-142",
      status: "success" as const,
    },
    {
      id: "seq2",
      prediction: "NifH Group III",
      confidence: "0.88",
      evalue: "2e-98",
      status: "success" as const,
    },
    {
      id: "seq3",
      prediction: "AnfH",
      confidence: "0.91",
      evalue: "4e-117",
      status: "success" as const,
    },
    {
      id: "seq4",
      prediction: "Non-nitrogenase",
      confidence: "—",
      evalue: "0.43",
      status: "error" as const,
    },
    {
      id: "seq5",
      prediction: "VnfH",
      confidence: "0.84",
      evalue: "6e-89",
      status: "success" as const,
    },
  ]

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
          <div style={styles.cardTitle}>Classification output</div>
          <button
            style={styles.btn}
            onClick={() => alert("Connect to your CSV download endpoint")}
          >
            Download CSV ↓
          </button>
        </div>

        <table style={styles.table}>
          <thead>
            <tr>
              {[
                "Sequence ID",
                "Prediction",
                "Confidence",
                "E-value",
                "Status",
              ].map((h) => (
                <th key={h} style={styles.th}>
                  {h}
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {mockResults.map((row) => (
              <tr
                key={row.id}
                style={{
                  background:
                    row.status === "error"
                      ? "rgba(255,92,92,0.03)"
                      : "transparent",
                }}
              >
                <td style={styles.td}>{row.id}</td>
                <td
                  style={{
                    ...styles.td,
                    color: row.status === "error" ? theme.red : theme.text,
                  }}
                >
                  {row.prediction}
                </td>
                <td style={styles.td}>{row.confidence}</td>
                <td style={styles.td}>{row.evalue}</td>
                <td style={styles.td}>
                  <span style={badgeStyle(row.status)}>
                    {row.status === "success"
                      ? "classified"
                      : "non-nitrogenase"}
                  </span>
                </td>
              </tr>
            ))}
          </tbody>
        </table>
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
          [
            "Non-protein sequences without Prodigal",
            "If submitting DNA/nucleotide sequences, enable the Prodigal CDS prediction option.",
          ],
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
  const [view, setView] = useState<View>("upload")
  const [jobData, setJobData] = useState<JobData | null>(null)

  const handleSubmit = (data: JobData) => {
    setJobData(data)
    setView("waiting")
  }

  return (
    <div style={styles.page}>
      <header style={styles.header}>
        <div style={styles.logo}>DiazoDB</div>
        <div style={styles.logoSub}>nitrogenase sequence classifier</div>
      </header>

      {view === "upload" && <UploadPage onSubmit={handleSubmit} />}
      {view === "waiting" && (
        <WaitingPage
          jobId={jobData?.jobId ?? "demo-job-id"}
          onSuccess={() => setView("results")}
          onFailure={() => setView("failure")}
        />
      )}
      {view === "results" && (
        <ResultsPage
          jobId={jobData?.jobId ?? ""}
          onReset={() => setView("upload")}
        />
      )}
      {view === "failure" && (
        <FailurePage
          jobId={jobData?.jobId ?? ""}
          onReset={() => setView("upload")}
        />
      )}
    </div>
  )
}
