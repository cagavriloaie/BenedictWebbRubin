#pragma codepage 65001
; =============================================================================
; BWRSflow Gas Flow Calculator - Inno Setup Script
; ELCOST Impex  |  v3.0 / 2026
;
; Pasi inainte de compilare:
;   1. Construiti configuratia Release x64 din Visual Studio (sau CMake):
;        cmake --preset x64-release
;        cmake --build ..\Out\Build\x64-release --config Release
;   2. Deschideti acest fisier in Inno Setup Compiler si apasati Build.
;   3. Installerul apare in folderul  Installer\  ca  BWRSflow_Setup_v3.0.exe
; =============================================================================

#define AppName      "BWRSflow Gas Flow Calculator"
#define AppVersion   "3.0"
#define AppPublisher "ELCOST Impex"
#define AppExeName   "BWRSflow.exe"

[Setup]
AppId={{A1B2C3D4-E5F6-7890-ABCD-EF1234567890}
AppName={#AppName}
AppVersion={#AppVersion}
AppVerName={#AppName} v{#AppVersion}
AppPublisher={#AppPublisher}
DefaultDirName={autopf}\{#AppPublisher}\BWRSflow
DefaultGroupName={#AppPublisher}
AllowNoIcons=yes
OutputDir=.
OutputBaseFilename=BWRSflow_Setup_v{#AppVersion}
Compression=lzma2/ultra64
SolidCompression=yes
WizardStyle=modern
ArchitecturesInstallIn64BitMode=x64
PrivilegesRequired=admin
SetupIconFile=..\App\BWRSflow.ico
UninstallDisplayName={#AppName}
UninstallDisplayIcon={app}\{#AppExeName}
MinVersion=10.0

[Languages]
Name: "romanian"; MessagesFile: "compiler:Languages\Romanian.isl"

[CustomMessages]
RunAfterInstall=Porneste {#AppName}

[Tasks]
Name: "desktopicon"; \
    Description: "Creaza pictograma pe &Desktop"; \
    GroupDescription: "{cm:AdditionalIcons}"; \
    Flags: unchecked

[Files]
Source: "..\Out\Build\x64-release\BWRSflow.exe"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\App\BWRSflow.ico";                    DestDir: "{app}"; Flags: ignoreversion

[Icons]
Name: "{group}\{#AppName}";    Filename: "{app}\{#AppExeName}"
Name: "{group}\{cm:UninstallProgram,{#AppName}}"; Filename: "{uninstallexe}"
Name: "{commondesktop}\{#AppName}"; Filename: "{app}\{#AppExeName}"; Tasks: desktopicon

[Run]
Filename: "{app}\{#AppExeName}"; \
    Description: "{cm:RunAfterInstall}"; \
    Flags: nowait postinstall skipifsilent
