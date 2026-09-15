import { execFileSync } from 'node:child_process'
import { readFileSync } from 'node:fs'
import { fileURLToPath } from 'node:url'
import { resolve } from 'node:path'

const frontendRoot = resolve(fileURLToPath(new URL('..', import.meta.url)))
const sourceRoot = resolve(frontendRoot, '..')
const expected = process.env.VERCEL_GIT_COMMIT_SHA?.trim()
  ?? execFileSync('git', ['rev-parse', 'HEAD'], { cwd: sourceRoot, encoding: 'utf8' }).trim()
const buildInfoPath = resolve(frontendRoot, 'dist', 'build-info.json')
const buildInfo = JSON.parse(readFileSync(buildInfoPath, 'utf8'))

if (buildInfo.sha !== expected) {
  throw new Error(`build-info SHA mismatch: expected ${expected}, got ${buildInfo.sha ?? '<missing>'}`)
}

if (!['git', 'vercel'].includes(buildInfo.source)) {
  throw new Error(`build-info source is invalid: ${buildInfo.source ?? '<missing>'}`)
}

console.log(`build-info verified: ${buildInfo.sha} (${buildInfo.source})`)
