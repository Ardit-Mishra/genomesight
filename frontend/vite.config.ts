import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'
import { execFileSync } from 'node:child_process'
import { fileURLToPath, URL } from 'node:url'

const sourceRoot = fileURLToPath(new URL('..', import.meta.url))

function sourceCommit() {
  const vercelCommit = process.env.VERCEL_GIT_COMMIT_SHA?.trim()
  if (vercelCommit) {
    return { sha: vercelCommit, source: 'vercel' }
  }

  const sha = execFileSync('git', ['rev-parse', 'HEAD'], {
    cwd: sourceRoot,
    encoding: 'utf8',
  }).trim()
  return { sha, source: 'git' }
}

function buildInfoPlugin() {
  const build = sourceCommit()
  return {
    name: 'genomesight-build-info',
    generateBundle() {
      this.emitFile({
        type: 'asset',
        fileName: 'build-info.json',
        source: `${JSON.stringify(build)}\n`,
      })
    },
  }
}

// https://vitejs.dev/config/
export default defineConfig({
  plugins: [react(), buildInfoPlugin()],
  server: {
    port: 3000,
    proxy: {
      '/api': {
        target: 'http://localhost:8000',
        changeOrigin: true,
      }
    }
  },
  resolve: {
    alias: { '@': fileURLToPath(new URL('./src', import.meta.url)) },
  }
})
