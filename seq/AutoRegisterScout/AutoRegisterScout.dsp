# Microsoft Developer Studio Project File - Name="AutoRegisterScout" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Dynamic-Link Library" 0x0102

CFG=AutoRegisterScout - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE
!MESSAGE NMAKE /f "AutoRegisterScout.mak".
!MESSAGE
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE
!MESSAGE NMAKE /f "AutoRegisterScout.mak" CFG="AutoRegisterScout - Win32 Debug"
!MESSAGE
!MESSAGE Possible choices for configuration are:
!MESSAGE
!MESSAGE "AutoRegisterScout - Win32 Debug" (based on "Win32 (x86) Dynamic-Link Library")
!MESSAGE

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName "AutoRegisterScout"
# PROP Scc_LocalPath ".."
CPP=xicl6.exe
MTL=midl.exe
RSC=rc.exe
# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 2
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "./Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MTd /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "AUTOREGISTERSCOUT_EXPORTS" /YX /FD /GZ /c
# ADD CPP /nologo /G6 /MDd /W3 /Gm /GR /GX /ZI /Od /D "GRE" /D "BUILD_SEQU" /D "MFCDebug" /D "_DEBUG" /D "RWDEBUG" /D "WIN32" /D "_MBCS" /D "_WINDOWS" /D "WinNT400" /D "MSVC60" /D AFX_NOVTABLE= /D "DEBUG" /D "_RWTOOLSDLL" /D "_UNICODE" /D "UNICODE" /D "O_DLL_PCLASS_ONLY" /D "_CONSOLE" /D "CSA_HAS_DLL" /D "ACE_HAS_DLL" /D "_WINDLL" /D "_AFXDLL" /D "SUPPORT_CT" /D "SUPPORT_PACE" /D "SUPPORT_iPAT" /D "NOCRYPT" /FR /FD /GZ /c
# SUBTRACT CPP /YX /Yc /Yu
# ADD BASE MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG" /d "_AFXDLL"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=xilink6.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /dll /debug /machine:I386 /pdbtype:sept
# ADD LINK32 rwtoold.lib aced.lib liboscfe.lib mscchk.lib mscvcch.lib i18nd.lib MrProtd.lib Sequenced.lib libMESd.lib MeasSectionsd.lib PerProxiesd.lib MeasPatientd.lib libSBBd.lib libSSLd.lib libRTd.lib libGSLd.lib SeqBufferd.lib libSeqUTd.lib libSeqUtild.lib MdhProxyd.lib libUICtrld.lib MrTraced.lib CoilIfd.lib MeasNucleiBased.lib libPACEd.lib libUICtrld.lib libUILinkd.lib /nologo /dll /debug /machine:I386 /out:"z:\n4\x86\prod\bin/AutoRegisterScoutd.dll" /pdbtype:sept
# Begin Target

# Name "AutoRegisterScout - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\AutoRegisterScout.cpp
# End Source File
# Begin Source File

SOURCE=..\common\CT\CT.cpp
# End Source File
# Begin Source File

SOURCE=..\common\CT\CT_UI.cpp
# End Source File
# Begin Source File

SOURCE=..\common\iPAT\iPAT.cpp
# End Source File
# Begin Source File

SOURCE=.\LocalSeqLoop.cpp
# End Source File
# Begin Source File

SOURCE=..\Kernels\SBBGREKernel.cpp
# End Source File
# Begin Source File

SOURCE=..\Kernels\SBBPhaseEncode.cpp
# End Source File
# Begin Source File

SOURCE=..\Kernels\SBBReadOut.cpp
# End Source File
# Begin Source File

SOURCE=..\..\AutoRegisterApply\AutoRegisterApply.cpp
# End Source File
# Begin Source File

SOURCE=..\..\AutoRegisterApply\mgh_isometry\MGH_AffineTransform.cpp
# End Source File
# Begin Source File

SOURCE=..\..\AutoRegisterApply\mgh_isometry\MGH_Quaternion.cpp
# End Source File

# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\common\CT\CT.h
# End Source File
# Begin Source File

SOURCE=..\common\CT\CT_UI.h
# End Source File
# Begin Source File

SOURCE=..\common\iPAT\iPAT.h
# End Source File
# Begin Source File

SOURCE=..\..\AutoRegisterApply\AutoRegisterApply.h
# End Source File
# Begin Source File

SOURCE=..\..\AutoRegisterApply\mgh_isometry\MGH_AffineTransform.h
# End Source File
# Begin Source File

SOURCE=..\..\AutoRegisterApply\mgh_isometry\MGH_Quaternion.h
# End Source File

# Begin Source File

SOURCE=.\LocalSeqLoop.h
# End Source File
# Begin Source File

SOURCE=..\Kernels\SBBGREKernel.h
# End Source File
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project
