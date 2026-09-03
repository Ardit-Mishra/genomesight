/**
 * Client-side result export.
 *
 * Exports run in the browser rather than through the API for three reasons: the
 * data is already here, it costs no round-trip, and it still works while a free-tier
 * backend is cold-starting -- which is exactly when a user is most likely to have
 * results on screen and least likely to tolerate another 50-second wait.
 */

/** RFC 4180: quote a field, and double any embedded quote. */
function csvCell(value: unknown): string {
  const text = value === null || value === undefined ? "" : String(value)
  return /[",\n\r]/.test(text) ? `"${text.replace(/"/g, '""')}"` : text
}

export function toCsv(headers: string[], rows: unknown[][]): string {
  return [headers, ...rows].map((row) => row.map(csvCell).join(",")).join("\r\n")
}

/**
 * Hand the file to the browser. Revoking the object URL matters: without it every
 * export leaks a blob for the lifetime of the tab.
 */
export function download(filename: string, content: string, mime = "text/csv;charset=utf-8"): void {
  const url = URL.createObjectURL(new Blob([content], { type: mime }))
  const anchor = document.createElement("a")
  anchor.href = url
  anchor.download = filename
  document.body.appendChild(anchor)
  anchor.click()
  anchor.remove()
  URL.revokeObjectURL(url)
}

/** Timestamped so repeated exports do not overwrite each other in the download folder. */
export function stamped(base: string, extension: string): string {
  const now = new Date().toISOString().replace(/[:.]/g, "-").slice(0, 19)
  return `genomesight-${base}-${now}.${extension}`
}
