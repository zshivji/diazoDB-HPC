import { expect, test } from "@playwright/test"
import { randomEmail } from "./utils/random"

// Classify page is public — no login required
test.use({ storageState: { cookies: [], origins: [] } })

const FASTA_SEQUENCE = `>seq1
MRQCIAIGYKGGIGKSTTSNTLAAMILAQNGAKVGLTGDAEGIDLGKIVEII
GHISSQYGKTIVNVGGPEEGLKFADLIQEGRIDKVEDLEKELETMLEKNLRA`

const FASTA_SEQUENCE_2 = `>seq2
ANRQCIAIGYKGGIGKSTTSNTLAAMILAQNGAKVGLTGDAEGIDLGKIVEI
MRQCIAIGYKGGIGKSTTSNTLAAMILAQNGAKVGLTGDAEGIDLGKIVEII`

// ---- Page load ----

test("Classify page is accessible and shows correct title", async ({
  page,
}) => {
  await page.goto("/classify")

  await expect(
    page.getByRole("heading", { name: "Nitrogenase Annotation and Classification" }),
  ).toBeVisible()
})

test("DiazoDB logo and subtitle are visible", async ({ page }) => {
  await page.goto("/classify")

  await expect(page.getByText("DiazoDB", { exact: true }).first()).toBeVisible()
  await expect(page.getByText("nitrogenase sequence classifier")).toBeVisible()
})

test("Citing banner is visible", async ({ page }) => {
  await page.goto("/classify")

  await expect(page.getByText("Please cite:")).toBeVisible()
})

test("Submit button is visible", async ({ page }) => {
  await page.goto("/classify")

  await expect(page.getByRole("button", { name: /Submit/i })).toBeVisible()
})

// ---- Form inputs ----

test("Email input is visible and editable", async ({ page }) => {
  await page.goto("/classify")

  const emailInput = page.getByPlaceholder("you@institution.edu")
  await expect(emailInput).toBeVisible()
  await expect(emailInput).toBeEditable()
})

test("Sequence textarea is visible and editable", async ({ page }) => {
  await page.goto("/classify")

  const textarea = page.getByPlaceholder(/>seq1/)
  await expect(textarea).toBeVisible()
  await expect(textarea).toBeEditable()
})

test("Prodigal checkbox is visible and unchecked by default", async ({
  page,
}) => {
  await page.goto("/classify")

  const checkbox = page.getByRole("checkbox")
  await expect(checkbox).toBeVisible()
  await expect(checkbox).not.toBeChecked()
})

test("File drop zone is visible", async ({ page }) => {
  await page.goto("/classify")

  await expect(
    page.getByText("Drop a FASTA file here, or click to browse"),
  ).toBeVisible()
})

// ---- Validation ----

test("Shows error when submitting without sequences", async ({ page }) => {
  await page.goto("/classify")

  await page.getByPlaceholder("you@institution.edu").fill(randomEmail())
  await page.getByRole("button", { name: /Submit/i }).click()

  await expect(
    page.getByText("Please provide sequences or upload a FASTA file"),
  ).toBeVisible()
})

test("Shows error when submitting without email", async ({ page }) => {
  await page.goto("/classify")

  await page.getByPlaceholder(/>seq1/).fill(FASTA_SEQUENCE)
  await page.getByRole("button", { name: /Submit/i }).click()

  await expect(
    page.getByText("Please provide an email address for results notification"),
  ).toBeVisible()
})

// ---- Manual paste submission ----

test.describe("Manual sequence paste", () => {
  test("Pasted sequences are accepted and show waiting page", async ({
    page,
  }) => {
    await page.goto("/classify")

    await page.getByPlaceholder("you@institution.edu").fill(randomEmail())
    await page.getByPlaceholder(/>seq1/).fill(FASTA_SEQUENCE)
    await page.getByRole("button", { name: /Submit/i }).click()

    await expect(page.getByText("Processing")).toBeVisible()
    const jobIdText = page.getByText(
      /[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}/i,
    )
    await expect(jobIdText).toBeVisible()
  })

  test("Multiple pasted sequences are accepted", async ({ page }) => {
    await page.goto("/classify")

    await page.getByPlaceholder("you@institution.edu").fill(randomEmail())
    await page
      .getByPlaceholder(/>seq1/)
      .fill(`${FASTA_SEQUENCE}\n\n${FASTA_SEQUENCE_2}`)
    await page.getByRole("button", { name: /Submit/i }).click()

    await expect(page.getByText("Processing")).toBeVisible()
  })

  test("Waiting page shows elapsed time after paste submission", async ({
    page,
  }) => {
    await page.goto("/classify")

    await page.getByPlaceholder("you@institution.edu").fill(randomEmail())
    await page.getByPlaceholder(/>seq1/).fill(FASTA_SEQUENCE)
    await page.getByRole("button", { name: /Submit/i }).click()

    await expect(page.getByText("Elapsed")).toBeVisible()
  })

  test("Simulate success after paste leads to results page", async ({
    page,
  }) => {
    await page.goto("/classify")

    await page.getByPlaceholder("you@institution.edu").fill(randomEmail())
    await page.getByPlaceholder(/>seq1/).fill(FASTA_SEQUENCE)
    await page.getByRole("button", { name: /Submit/i }).click()

    await expect(page.getByText("Processing")).toBeVisible()
    await page.getByRole("button", { name: /Simulate success/i }).click()

    await expect(page.getByRole("heading", { name: "Results" })).toBeVisible()
    await expect(page.getByText("Classification output")).toBeVisible()
    await expect(
      page.getByRole("button", { name: /Download CSV/i }),
    ).toBeVisible()
  })

  test("Simulate failure after paste leads to failure page", async ({
    page,
  }) => {
    await page.goto("/classify")

    await page.getByPlaceholder("you@institution.edu").fill(randomEmail())
    await page.getByPlaceholder(/>seq1/).fill(FASTA_SEQUENCE)
    await page.getByRole("button", { name: /Submit/i }).click()

    await expect(page.getByText("Processing")).toBeVisible()
    await page.getByRole("button", { name: /Simulate failure/i }).click()

    await expect(page.getByText("Job failed")).toBeVisible()
    await expect(page.getByText("Common causes")).toBeVisible()
  })
})

// ---- File upload submission ----

test.describe("File upload", () => {
  test("Uploading a FASTA file shows filename in drop zone", async ({
    page,
  }) => {
    await page.goto("/classify")

    const fileInput = page.locator('input[type="file"]')
    await fileInput.setInputFiles({
      name: "sequences.fasta",
      mimeType: "text/plain",
      buffer: Buffer.from(FASTA_SEQUENCE),
    })

    await expect(page.getByText("✓ sequences.fasta")).toBeVisible()
  })

  test("File upload submission shows waiting page with job ID", async ({
    page,
  }) => {
    await page.goto("/classify")

    await page.getByPlaceholder("you@institution.edu").fill(randomEmail())

    const fileInput = page.locator('input[type="file"]')
    await fileInput.setInputFiles({
      name: "sequences.fasta",
      mimeType: "text/plain",
      buffer: Buffer.from(FASTA_SEQUENCE),
    })

    await expect(page.getByText("✓ sequences.fasta")).toBeVisible()
    await page.getByRole("button", { name: /Submit/i }).click()

    await expect(page.getByText("Processing")).toBeVisible()
    const jobIdText = page.getByText(
      /[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}/i,
    )
    await expect(jobIdText).toBeVisible()
  })

  test("Multi-sequence FASTA file upload is accepted", async ({ page }) => {
    await page.goto("/classify")

    await page.getByPlaceholder("you@institution.edu").fill(randomEmail())

    const fileInput = page.locator('input[type="file"]')
    await fileInput.setInputFiles({
      name: "multi.fasta",
      mimeType: "text/plain",
      buffer: Buffer.from(`${FASTA_SEQUENCE}\n\n${FASTA_SEQUENCE_2}`),
    })

    await expect(page.getByText("✓ multi.fasta")).toBeVisible()
    await page.getByRole("button", { name: /Submit/i }).click()

    await expect(page.getByText("Processing")).toBeVisible()
  })

  test("Simulate success after file upload leads to results page", async ({
    page,
  }) => {
    await page.goto("/classify")

    await page.getByPlaceholder("you@institution.edu").fill(randomEmail())

    const fileInput = page.locator('input[type="file"]')
    await fileInput.setInputFiles({
      name: "sequences.fasta",
      mimeType: "text/plain",
      buffer: Buffer.from(FASTA_SEQUENCE),
    })

    await page.getByRole("button", { name: /Submit/i }).click()
    await expect(page.getByText("Processing")).toBeVisible()
    await page.getByRole("button", { name: /Simulate success/i }).click()

    await expect(page.getByRole("heading", { name: "Results" })).toBeVisible()
    await expect(page.getByText("Classification output")).toBeVisible()
  })

  test("Simulate failure after file upload leads to failure page", async ({
    page,
  }) => {
    await page.goto("/classify")

    await page.getByPlaceholder("you@institution.edu").fill(randomEmail())

    const fileInput = page.locator('input[type="file"]')
    await fileInput.setInputFiles({
      name: "sequences.fasta",
      mimeType: "text/plain",
      buffer: Buffer.from(FASTA_SEQUENCE),
    })

    await page.getByRole("button", { name: /Submit/i }).click()
    await expect(page.getByText("Processing")).toBeVisible()
    await page.getByRole("button", { name: /Simulate failure/i }).click()

    await expect(page.getByText("Job failed")).toBeVisible()
  })

  test("File upload without email shows validation error", async ({ page }) => {
    await page.goto("/classify")

    const fileInput = page.locator('input[type="file"]')
    await fileInput.setInputFiles({
      name: "sequences.fasta",
      mimeType: "text/plain",
      buffer: Buffer.from(FASTA_SEQUENCE),
    })

    await page.getByRole("button", { name: /Submit/i }).click()

    await expect(
      page.getByText(
        "Please provide an email address for results notification",
      ),
    ).toBeVisible()
  })

  test("Uploading a file populates textarea content", async ({ page }) => {
    await page.goto("/classify")

    const fileInput = page.locator('input[type="file"]')
    await fileInput.setInputFiles({
      name: "sequences.fasta",
      mimeType: "text/plain",
      buffer: Buffer.from(FASTA_SEQUENCE),
    })

    // After upload, the textarea should contain the file contents
    const textarea = page.getByPlaceholder(/>seq1/)
    await expect(textarea).not.toBeEmpty()
  })
})

// ---- Navigation flows ----

test("Try again button on failure page returns to upload form", async ({
  page,
}) => {
  await page.goto("/classify")

  await page.getByPlaceholder("you@institution.edu").fill(randomEmail())
  await page.getByPlaceholder(/>seq1/).fill(FASTA_SEQUENCE)
  await page.getByRole("button", { name: /Submit/i }).click()

  await page.getByRole("button", { name: /Simulate failure/i }).click()
  await page.getByRole("button", { name: /Try again/i }).click()

  await expect(
    page.getByRole("heading", { name: "Classify Nitrogenase Sequences" }),
  ).toBeVisible()
})

test("Submit another job button on results page returns to upload form", async ({
  page,
}) => {
  await page.goto("/classify")

  await page.getByPlaceholder("you@institution.edu").fill(randomEmail())
  await page.getByPlaceholder(/>seq1/).fill(FASTA_SEQUENCE)
  await page.getByRole("button", { name: /Submit/i }).click()

  await page.getByRole("button", { name: /Simulate success/i }).click()
  await page.getByRole("button", { name: /Submit another job/i }).click()

  await expect(
    page.getByRole("heading", { name: "Classify Nitrogenase Sequences" }),
  ).toBeVisible()
})

// ---- Prodigal option ----

test("Prodigal checkbox can be toggled", async ({ page }) => {
  await page.goto("/classify")

  const checkbox = page.getByRole("checkbox")
  await expect(checkbox).not.toBeChecked()

  await page
    .locator("label")
    .filter({ hasText: "Predict CDS with Prodigal" })
    .click()
  await expect(checkbox).toBeChecked()

  await page.getByText("Predict CDS with Prodigal").click()
  await expect(checkbox).not.toBeChecked()
})
