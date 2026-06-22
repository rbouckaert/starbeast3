#!/bin/bash
#
# Build and release a BEAST 2.8 package.
#
# Runs `mvn package` to produce the BEAST package ZIP (via maven-assembly-plugin),
# then optionally creates a GitHub release with the ZIP attached.
#
# Usage:
#   ./release.sh              # build + create package ZIP
#   ./release.sh --release    # also create a GitHub release with the ZIP attached
#
# The ZIP can then be submitted to CBAN (CompEvol/CBAN) by adding an entry
# to packages2.8.xml via pull request. See README.md for details.
#
set -euo pipefail

# --- Extract metadata from version.xml ---

if [[ ! -f version.xml ]]; then
    echo "ERROR: version.xml not found in $(pwd)" >&2
    exit 1
fi

# Parse only the <package> element (not <provider classname="..."> etc.)
PKG_LINE=$(grep '<package ' version.xml | head -1)
PKG_NAME=$(echo "$PKG_LINE" | sed "s/.*name=['\"]\\([^'\"]*\\)['\"].*/\\1/")
VERSION=$(echo "$PKG_LINE" | sed "s/.*version=['\"]\\([^'\"]*\\)['\"].*/\\1/")

if [[ -z "$PKG_NAME" || -z "$VERSION" ]]; then
    echo "ERROR: could not parse name/version from version.xml" >&2
    exit 1
fi

ZIP_NAME="${PKG_NAME}.v${VERSION}.zip"
GITHUB_REPO=$(git remote get-url origin 2>/dev/null \
    | sed 's|.*github.com[:/]\(.*\)\.git$|\1|; s|.*github.com[:/]\(.*\)$|\1|')

echo "=== ${PKG_NAME} v${VERSION} ==="
echo ""

# --- Step 1: Maven build + assemble package ZIP ---

echo "--- Building with Maven ---"
mvn clean package -DskipTests
echo ""

# Locate the ZIP produced by maven-assembly-plugin
ZIP_PATH="target/${ZIP_NAME}"
if [[ ! -f "$ZIP_PATH" ]]; then
    echo "ERROR: expected ZIP not found at ${ZIP_PATH}" >&2
    exit 1
fi

# Copy to project root for convenience
cp "$ZIP_PATH" "$ZIP_NAME"

echo "=== Package: ${ZIP_NAME} ==="
unzip -l "$ZIP_NAME"

# --- Step 2: Optionally create GitHub release ---

if [[ "${1:-}" == "--release" ]]; then
    if [[ -z "$GITHUB_REPO" ]]; then
        echo "ERROR: could not determine GitHub repo from git remote" >&2
        exit 1
    fi

    echo ""
    echo "--- Creating GitHub release v${VERSION} on ${GITHUB_REPO} ---"
    gh release create "v${VERSION}" "$ZIP_NAME" \
        --repo "$GITHUB_REPO" \
        --title "${PKG_NAME} v${VERSION}" \
        --generate-notes

    DOWNLOAD_URL="https://github.com/${GITHUB_REPO}/releases/download/v${VERSION}/${ZIP_NAME}"
    echo ""
    echo "=== Release created ==="
    echo "URL: https://github.com/${GITHUB_REPO}/releases/tag/v${VERSION}"
    echo ""
    echo "--- Next step: submit to CBAN ---"
    echo "Add this entry to packages2.8.xml in https://github.com/CompEvol/CBAN via pull request:"
    echo ""
    cat <<XMLEOF
    <package name="${PKG_NAME}" version="${VERSION}"
        url="${DOWNLOAD_URL}"
        projectURL="https://github.com/${GITHUB_REPO}"
        description="Optimised relaxed clock model">
        <depends on="BEAST.base" atleast="2.8.0"/>
    </package>
XMLEOF
fi