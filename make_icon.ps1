# Genereaza BWR.ico — contor industrial (manometru cu ac)
# Ruleaza o singura data din directorul proiectului:
#   powershell -ExecutionPolicy Bypass -File make_icon.ps1

Add-Type -AssemblyName System.Drawing

function New-GaugeIcon {
    param([int]$Size)

    $bmp = New-Object System.Drawing.Bitmap($Size, $Size,
        [System.Drawing.Imaging.PixelFormat]::Format32bppArgb)
    $g = [System.Drawing.Graphics]::FromImage($bmp)
    $g.SmoothingMode     = [System.Drawing.Drawing2D.SmoothingMode]::AntiAlias
    $g.TextRenderingHint = [System.Drawing.Text.TextRenderingHint]::AntiAliasGridFit

    $cx = [float]($Size / 2)
    $cy = [float]($Size / 2)

    # ── fundal ───────────────────────────────────────────────────────────────
    $bg = New-Object System.Drawing.SolidBrush(
        [System.Drawing.Color]::FromArgb(255, 8, 25, 62))
    $g.FillEllipse($bg, 0, 0, $Size, $Size)
    $bg.Dispose()

    # strat interior luminos
    $glow = New-Object System.Drawing.SolidBrush(
        [System.Drawing.Color]::FromArgb(60, 45, 110, 200))
    $gs = [float]($Size * 0.68)
    $go = [float](($Size - $gs) / 2)
    $g.FillEllipse($glow, $go, $go, $gs, $gs)
    $glow.Dispose()

    # ── inel exterior ────────────────────────────────────────────────────────
    $rw = [Math]::Max(1.0, $Size * 0.048)
    $rp = New-Object System.Drawing.Pen(
        [System.Drawing.Color]::FromArgb(255, 85, 145, 215), $rw)
    $m = [float]($rw / 2)
    $g.DrawEllipse($rp, $m, $m, [float]($Size - $rw), [float]($Size - $rw))
    $rp.Dispose()

    # ── inel interior (cadran) ───────────────────────────────────────────────
    $rw2 = [Math]::Max(0.5, $Size * 0.02)
    $rp2 = New-Object System.Drawing.Pen(
        [System.Drawing.Color]::FromArgb(120, 80, 130, 200), $rw2)
    $ir  = [float]($cx * 0.88)
    $g.DrawEllipse($rp2, [float]($cx - $ir), [float]($cy - $ir),
        [float]($ir * 2), [float]($ir * 2))
    $rp2.Dispose()

    # ── arc scala (225° — 315° sens orar, deschis jos) ───────────────────────
    if ($Size -ge 32) {
        $arcW = [Math]::Max(1.0, $Size * 0.025)
        $ap   = New-Object System.Drawing.Pen(
            [System.Drawing.Color]::FromArgb(180, 100, 160, 230), $arcW)
        $ar   = [float]($cx * 0.76)
        $g.DrawArc($ap,
            [float]($cx - $ar), [float]($cy - $ar),
            [float]($ar * 2), [float]($ar * 2),
            135.0, 270.0)
        $ap.Dispose()
    }

    # ── tick marks (12 pozitii, din 30 in 30 grade) ──────────────────────────
    $rOuter = [float]($cx * 0.82)
    for ($i = 0; $i -lt 12; $i++) {
        $rad    = [float]($i * 30.0 * [Math]::PI / 180.0 - [Math]::PI / 2.0)
        $major  = ($i % 3 -eq 0)
        $rInner = if ($major) { [float]($cx * 0.62) } else { [float]($cx * 0.72) }
        $tw     = if ($major) {
            [Math]::Max(1.0, $Size * 0.032)
        } else {
            [Math]::Max(0.5, $Size * 0.016)
        }
        $col = if ($major) {
            [System.Drawing.Color]::FromArgb(255, 215, 228, 255)
        } else {
            [System.Drawing.Color]::FromArgb(170, 150, 175, 215)
        }
        $tp = New-Object System.Drawing.Pen($col, $tw)
        $tp.StartCap = [System.Drawing.Drawing2D.LineCap]::Round
        $tp.EndCap   = [System.Drawing.Drawing2D.LineCap]::Round
        $g.DrawLine($tp,
            [float]($cx + $rOuter * [Math]::Cos($rad)),
            [float]($cy + $rOuter * [Math]::Sin($rad)),
            [float]($cx + $rInner * [Math]::Cos($rad)),
            [float]($cy + $rInner * [Math]::Sin($rad)))
        $tp.Dispose()
    }

    # ── ac indicator (la ora 2, aprox 2/3 din scala) ─────────────────────────
    $nRad  = [float](60.0 * [Math]::PI / 180.0 - [Math]::PI / 2.0)
    $nLen  = [float]($cx * 0.57)
    $tLen  = [float]($cx * 0.17)
    $nw    = [Math]::Max(1.5, $Size * 0.042)
    $np    = New-Object System.Drawing.Pen(
        [System.Drawing.Color]::FromArgb(255, 255, 72, 50), $nw)
    $np.StartCap = [System.Drawing.Drawing2D.LineCap]::Round
    $np.EndCap   = [System.Drawing.Drawing2D.LineCap]::Round
    $g.DrawLine($np,
        [float]($cx - $tLen * [Math]::Cos($nRad)),
        [float]($cy - $tLen * [Math]::Sin($nRad)),
        [float]($cx + $nLen * [Math]::Cos($nRad)),
        [float]($cy + $nLen * [Math]::Sin($nRad)))
    $np.Dispose()

    # ── buton central ─────────────────────────────────────────────────────────
    $hr = [float]([Math]::Max(2.0, $cx * 0.12))
    $hb = New-Object System.Drawing.SolidBrush(
        [System.Drawing.Color]::FromArgb(255, 185, 205, 240))
    $g.FillEllipse($hb,
        [float]($cx - $hr), [float]($cy - $hr),
        [float]($hr * 2), [float]($hr * 2))
    $hb.Dispose()
    $hbo = New-Object System.Drawing.SolidBrush(
        [System.Drawing.Color]::FromArgb(120, 20, 60, 130))
    $hro = [float]($hr * 0.55)
    $g.FillEllipse($hbo,
        [float]($cx - $hro), [float]($cy - $hro),
        [float]($hro * 2), [float]($hro * 2))
    $hbo.Dispose()

    # ── text "BWR" (doar pentru dimensiuni >= 48px) ───────────────────────────
    if ($Size -ge 48) {
        $fs   = [float]($Size * 0.135)
        $font = New-Object System.Drawing.Font(
            'Arial', $fs,
            [System.Drawing.FontStyle]::Bold,
            [System.Drawing.GraphicsUnit]::Pixel)
        $tb   = New-Object System.Drawing.SolidBrush(
            [System.Drawing.Color]::FromArgb(200, 200, 225, 255))
        $tsz  = $g.MeasureString('BWR', $font)
        $g.DrawString('BWR', $font, $tb,
            [float]($cx - $tsz.Width / 2),
            [float]($cy + $cx * 0.30))
        $font.Dispose()
        $tb.Dispose()
    }

    $g.Dispose()
    return $bmp
}

function Save-IcoFile {
    param(
        [System.Drawing.Bitmap[]]$Bitmaps,
        [string]$Path
    )

    $pngList = [System.Collections.Generic.List[byte[]]]::new()
    foreach ($b in $Bitmaps) {
        $ms = New-Object System.IO.MemoryStream
        $b.Save($ms, [System.Drawing.Imaging.ImageFormat]::Png)
        $pngList.Add($ms.ToArray())
        $ms.Dispose()
    }

    $n      = $Bitmaps.Count
    $fs     = [System.IO.File]::OpenWrite($Path)
    $writer = New-Object System.IO.BinaryWriter($fs)

    # ICO header (6 bytes)
    $writer.Write([uint16]0)   # reserved
    $writer.Write([uint16]1)   # type = ICO
    $writer.Write([uint16]$n)  # numar imagini

    # offset prima imagine = header + n * 16 bytes director
    $offset = [uint32](6 + 16 * $n)

    # Director imagini (16 bytes / imagine)
    for ($i = 0; $i -lt $n; $i++) {
        $s = $Bitmaps[$i].Width
        $writer.Write([byte]$(if ($s -ge 256) { 0 } else { $s }))  # width
        $writer.Write([byte]$(if ($s -ge 256) { 0 } else { $s }))  # height
        $writer.Write([byte]0)      # color count
        $writer.Write([byte]0)      # rezervat
        $writer.Write([uint16]1)    # plane-uri culoare
        $writer.Write([uint16]32)   # biti/pixel
        $writer.Write([uint32]$pngList[$i].Length)  # dimensiune date
        $writer.Write([uint32]$offset)              # offset date
        $offset += [uint32]$pngList[$i].Length
    }

    # Date imagini (PNG embed)
    foreach ($data in $pngList) { $writer.Write($data) }

    $writer.Close()
    $fs.Dispose()
}

# Genereaza 4 dimensiuni standard pentru .ico
$sizes   = 16, 32, 48, 256
$bitmaps = @($sizes | ForEach-Object { New-GaugeIcon $_ })

$out = Join-Path $PSScriptRoot 'BWR.ico'
Save-IcoFile -Bitmaps $bitmaps -Path $out
foreach ($b in $bitmaps) { $b.Dispose() }

Write-Host "Creat: $out" -ForegroundColor Green
