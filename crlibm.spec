Summary: cr-libm, a portable, efficient, correctly rounded mathematical library.
Name: crlibm
Version: 0.8
Release: 1
License: GNU LESSER GENERAL PUBLIC LICENSE
Group: System/Mathematics
URL: http://lipforge.ens-lyon.fr/projects/crlibm
Source0: %{name}-%{version}.tar.gz
BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-buildroot

%description The goal of the crlibm project is to design a
mathematical library which offers proven correct rounding, with
performance and memory usage comparable to current libraries.

%prep
%setup -q

%build
%configure
make

%install
rm -rf $RPM_BUILD_ROOT
%makeinstall

%clean
rm -rf $RPM_BUILD_ROOT

%files
%defattr(-,root,root,-)
%{_bindir}/*
%{_libdir}/*
%{_includedir}/*




%changelog
* Thu Oct 28 2004 Philippe Defert <defert@lxcert-i386.cern.ch> 
- Initial build.


