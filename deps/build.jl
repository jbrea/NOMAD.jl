using BinDeps
using Compat

@BinDeps.setup

#TODO automatically detect version
nomadversion = "3.8.1"

libnomad = library_dependency("libnomad")
libsgtelib = library_dependency("libsgtelib")

filename = "unix_linux/NOMAD.zip"
if is_apple()
    filename = "mac/NOMAD.dmg"
elseif is_windows()
    filename = "windows/NOMAD_setup.exe"
end
provides(Sources,
         URI("https://sourceforge.net/projects/nomad-bb-opt/files/$filename"),
         libnomad, unpacked_dir = "nomad.$nomadversion")

depsdir = BinDeps.depsdir(libnomad)
installdir = joinpath(depsdir, "src", "nomad", "install")
libdir = joinpath(depsdir, "src", "nomad", "lib")
provides(BuildProcess,
         (@build_steps begin
            GetSources(libnomad)
            `mv src/nomad.$nomadversion src/nomad`
            @build_steps begin
                ChangeDirectory(installdir)
                `bash install.sh`
            end
        end), [libnomad, libsgtelib], os = :Linux, installed_libpath=libdir)

#TODO build for other OS

@BinDeps.install Dict(:libnomad => :libnomad, :libsgtelib => :libsgtelib)

