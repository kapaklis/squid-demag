#!/bin/sh
# Build "SQUID Demag.app" next to this script: a launcher that runs gui.py with the python3 found now.
# Drag the app to /Applications or the Dock. Rebuild if you move this folder or change Python.
set -e
cd "$(dirname "$0")"
DIR="$PWD"
PY="$(command -v python3)"
APP="SQUID Demag.app"

rm -rf "$APP"
mkdir -p "$APP/Contents/MacOS" "$APP/Contents/Resources"
cp icon.icns "$APP/Contents/Resources/icon.icns"

cat > "$APP/Contents/MacOS/squid-demag" <<EOF
#!/bin/sh
exec "$PY" "$DIR/gui.py" "\$@"
EOF
chmod +x "$APP/Contents/MacOS/squid-demag"

cat > "$APP/Contents/Info.plist" <<EOF
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
  <key>CFBundleName</key><string>SQUID Demag</string>
  <key>CFBundleDisplayName</key><string>SQUID Demag</string>
  <key>CFBundleIdentifier</key><string>se.uu.physics.squid-demag</string>
  <key>CFBundleExecutable</key><string>squid-demag</string>
  <key>CFBundleIconFile</key><string>icon</string>
  <key>CFBundlePackageType</key><string>APPL</string>
  <key>CFBundleShortVersionString</key><string>1.0</string>
  <key>NSHighResolutionCapable</key><true/>
</dict>
</plist>
EOF
echo "built $DIR/$APP (python: $PY)"
