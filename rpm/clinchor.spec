Summary:	Accelerator code
Name:		clinchor
License:	EPICS Open license http://www.aps.anl.gov/epics/license/open.php
Group:		Applications/Databases
URL:		http://www.aps.anl.gov/asd/oag/oaghome.shtml
Packager:	Robert Soliday <soliday@aps.anl.gov>
Prefix:		%{_bindir}
Autoreq:	0
Version:	2.0
Release:	1
Source:		clinchor-2.0.tar.gz


%define debug_package %{nil}
%undefine __check_files
%description
Binary package for Clinchor. Clinchor calculates the growth rates 
of longitudinal and transverse coupled bunch modes in an electron 
storage ring.

%prep
%setup

%build
%install
mkdir -p %{buildroot}%{_bindir}
install -s -m 755 clinchor %{buildroot}%{_bindir}/clinchor

%files

%{_bindir}/clinchor

