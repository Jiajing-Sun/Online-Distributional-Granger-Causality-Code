# GitHub Push Commands

The local repository is ready. To create a private GitHub repository and push it, authenticate GitHub CLI first:

```bash
gh auth login
```

Then run from this repository root:

```bash
gh repo create Online-Distributional-Granger-Causality-Code --private --source . --remote origin --push
```

If the repository already exists on GitHub, add the remote and push:

```bash
git remote add origin git@github.com:<YOUR-GITHUB-USERNAME>/Online-Distributional-Granger-Causality-Code.git
git push -u origin main
```

The repository is intentionally private/code-only. Do not add raw Deribit transaction data or generated output folders.
