#!/usr/bin/env bash
# =============================================================================
#  BWRSflow — install.sh
#  Descarca, compileaza si porneste BWRSflow pe Linux.
#  Utilizare:  chmod +x install.sh && ./install.sh
# =============================================================================

set -euo pipefail

REPO_URL="https://github.com/cagavriloaie/BenedictWebbRubin"
REPO_DIR="BenedictWebbRubin"
BUILD_DIR="Build"
BINARY="BWRSflow"

# ── Culori ────────────────────────────────────────────────────────────────────
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
CYAN='\033[0;36m'
BOLD='\033[0;37m'
RESET='\033[0m'

info()    { echo -e "${CYAN}[INFO]${RESET}  $*"; }
ok()      { echo -e "${GREEN}[OK]${RESET}    $*"; }
warn()    { echo -e "${YELLOW}[WARN]${RESET}  $*"; }
error()   { echo -e "${RED}[ERROR]${RESET} $*" >&2; exit 1; }
section() { echo -e "\n${BOLD}── $* ──────────────────────────────────────────${RESET}"; }

# ── Detectare distributie ─────────────────────────────────────────────────────
detect_distro() {
    if   command -v apt-get &>/dev/null; then echo "debian"
    elif command -v dnf     &>/dev/null; then echo "fedora"
    elif command -v pacman  &>/dev/null; then echo "arch"
    elif command -v zypper  &>/dev/null; then echo "opensuse"
    else echo "unknown"
    fi
}

# ── Instalare dependente ──────────────────────────────────────────────────────
install_deps() {
    section "Instalare dependente"
    local distro
    distro=$(detect_distro)

    case "$distro" in
        debian)
            info "Detectat: Ubuntu / Debian / Linux Mint"
            info "Rulare: sudo apt-get update && sudo apt-get install ..."
            sudo apt-get update -qq
            sudo apt-get install -y git build-essential cmake
            ;;
        fedora)
            info "Detectat: Fedora / RHEL / CentOS Stream"
            sudo dnf install -y git gcc-c++ cmake
            ;;
        arch)
            info "Detectat: Arch Linux / Manjaro"
            sudo pacman -S --noconfirm git gcc cmake
            ;;
        opensuse)
            info "Detectat: openSUSE"
            sudo zypper install -y git gcc-c++ cmake
            ;;
        *)
            warn "Distributie necunoscuta. Asigurati-va ca aveti instalat:"
            warn "  git, g++ (>= 10), cmake (>= 3.10)"
            ;;
    esac
}

# ── Verificare versiuni ───────────────────────────────────────────────────────
check_versions() {
    section "Verificare versiuni"

    # Git
    if command -v git &>/dev/null; then
        ok "git: $(git --version)"
    else
        error "git nu este instalat."
    fi

    # C++ compiler
    CXX_CMD=""
    if command -v g++ &>/dev/null; then
        CXX_CMD="g++"
    elif command -v clang++ &>/dev/null; then
        CXX_CMD="clang++"
    else
        error "Nu s-a gasit g++ sau clang++. Instalati build-essential sau gcc-c++."
    fi

    CXX_VER=$($CXX_CMD --version | head -1)
    ok "compilator: $CXX_VER"

    # Versiune minima GCC 10 / Clang 12
    if command -v g++ &>/dev/null; then
        GCC_MAJOR=$(g++ -dumpversion | cut -d. -f1)
        if [ "$GCC_MAJOR" -lt 10 ]; then
            error "GCC $GCC_MAJOR detectat — necesita >= 10 pentru C++20. Actualizati compilatorul."
        fi
    fi

    # CMake
    if command -v cmake &>/dev/null; then
        ok "cmake: $(cmake --version | head -1)"
    else
        error "cmake nu este instalat."
    fi
}

# ── Clonare / actualizare repo ────────────────────────────────────────────────
clone_repo() {
    section "Descarcare sursa"

    if [ -d "$REPO_DIR/.git" ]; then
        info "Repo existent gasit — actualizare cu git pull..."
        git -C "$REPO_DIR" pull
    else
        info "Clonare din $REPO_URL ..."
        git clone "$REPO_URL" "$REPO_DIR"
    fi

    ok "Sursa descarcata in: $(realpath "$REPO_DIR")"
}

# ── Compilare ─────────────────────────────────────────────────────────────────
build_project() {
    section "Compilare"

    cd "$REPO_DIR"

    info "Configurare CMake..."
    cmake -B "$BUILD_DIR" -S . -DCMAKE_BUILD_TYPE=Release

    info "Compilare..."
    cmake --build "$BUILD_DIR" --config Release

    if [ -f "$BUILD_DIR/$BINARY" ]; then
        ok "Executabil creat: $(realpath "$BUILD_DIR/$BINARY")"
    else
        error "Compilarea a esuat — executabilul nu a fost creat."
    fi
}

# ── Rulare ────────────────────────────────────────────────────────────────────
run_app() {
    section "Pornire aplicatie"
    info "Lansare $BINARY..."
    echo ""
    cd "$BUILD_DIR"
    ./$BINARY
}

# ── Main ──────────────────────────────────────────────────────────────────────
main() {
    echo -e "\n${YELLOW}══════════════════════════════════════════════════════${RESET}"
    echo -e "${YELLOW}  BWRSflow — Gas Flow Calculator  |  ELCOST Impex     ${RESET}"
    echo -e "${YELLOW}  Script instalare Linux  v3.0/2026                   ${RESET}"
    echo -e "${YELLOW}══════════════════════════════════════════════════════${RESET}\n"

    install_deps
    check_versions
    clone_repo
    build_project
    run_app
}

main "$@"
