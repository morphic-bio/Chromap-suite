#!/usr/bin/env python3
"""Build tested Chromap release artifacts. Publication is owned by CI."""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tarfile
import tempfile

ROOT = Path(__file__).resolve().parents[2]
BINARIES = ("chromap", "rapidmacs", "chromap_callpeaks", "chromap_lib_runner",
            "chromap_atac_spill_materializer")


def run(args, **kwargs):
    return subprocess.run([str(a) for a in args], check=True, **kwargs)


def output(args, **kwargs):
    return subprocess.check_output([str(a) for a in args], text=True, **kwargs).strip()


def suite_version():
    return re.search(r'#define CHROMAP_SUITE_VERSION "([^"]+)"',
                     (ROOT / "src/version.h").read_text())[1]


def version_info(value):
    value = value.removeprefix("v")
    match = re.fullmatch(r"(\d+\.\d+\.\d+)(?:-([A-Za-z0-9][A-Za-z0-9.-]*))?", value)
    if not match:
        raise ValueError(f"invalid release version: {value}")
    if match[1] != suite_version():
        raise ValueError(f"release {value} disagrees with src/version.h ({suite_version()})")
    return "v" + value, match[1], value.replace("-", "~", 1)


def provenance():
    snapshot = ROOT / ".release-source.json"
    if snapshot.exists():
        return json.loads(snapshot.read_text())
    return {
        "source_revision": output(["git", "rev-parse", "HEAD"], cwd=ROOT),
        "rapidmacs_revision": output(["git", "rev-parse", "HEAD:third_party/rapidmacs"], cwd=ROOT),
        "source_date_epoch": int(output(["git", "show", "-s", "--format=%ct", "HEAD"], cwd=ROOT)),
    }


def distribution():
    fields = dict(line.split("=", 1) for line in Path("/etc/os-release").read_text().splitlines()
                  if "=" in line)
    name, version = fields["ID"].strip('"'), fields["VERSION_ID"].strip('"')
    if name != "ubuntu" or version not in ("22.04", "24.04"):
        raise ValueError("release builds require Ubuntu 22.04 or 24.04; use the release build container")
    return name + version


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def runtime_dependencies(prefix):
    # dpkg's shlibs database supplies dependencies for this build baseline.
    with tempfile.TemporaryDirectory() as temp:
        temp = Path(temp)
        (temp / "debian").mkdir()
        shutil.copy2(ROOT / "debian/control", temp / "debian/control")
        text = output(["dpkg-shlibdeps", "-O", *[f"-e{prefix / 'bin' / b}" for b in BINARIES]], cwd=temp)
    return next(line.split("=", 1)[1] for line in text.splitlines()
                if line.startswith("shlibs:Depends="))


def stage_release(dest, version):
    tag, suite, deb_version = version_info(version)
    arch = output(["dpkg", "--print-architecture"])
    if arch != "amd64":
        raise ValueError("this release matrix currently supports amd64")
    distro = distribution()
    notes = ROOT / "docs" / f"RELEASE_NOTES_v{suite}.md"
    if not notes.is_file() or not notes.stat().st_size:
        raise ValueError(f"missing release notes: {notes}")
    got = output([ROOT / "chromap", "--version"], stderr=subprocess.STDOUT)
    if got != suite:
        raise ValueError(f"built binary version {got!r} != {suite!r}")
    dest.mkdir(parents=True, exist_ok=True)
    if any(dest.iterdir()):
        raise ValueError(f"staging directory must be empty: {dest}")
    for subdir in ("bin", "lib", "include/chromap-suite", "include/rapidmacs", "include/htslib",
                   "share/chromap-suite", "share/doc/chromap-suite", "share/licenses/chromap-suite"):
        (dest / subdir).mkdir(parents=True, exist_ok=True)
    for name in BINARIES:
        src = ROOT / name
        if not src.is_file() or not os.access(src, os.X_OK):
            raise ValueError(f"missing built executable: {src}")
        shutil.copy2(src, dest / "bin" / name)
    for src in (ROOT / "libchromap.a", ROOT / "third_party/rapidmacs/lib/librapidmacs.a"):
        if not src.is_file() or not src.stat().st_size:
            raise ValueError(f"missing built archive: {src}")
        shutil.copy2(src, dest / "lib" / src.name)
    for src in (ROOT / "src").rglob("*"):
        if src.is_file() and src.suffix in (".h", ".hpp"):
            target = dest / "include/chromap-suite" / src.relative_to(ROOT / "src")
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, target)
    for name, src in (("htslib", ROOT / "third_party/htslib/htslib"),
                      ("rapidmacs", ROOT / "third_party/rapidmacs/include/rapidmacs")):
        shutil.copytree(src, dest / "include" / name, dirs_exist_ok=True)
    for src, name in ((ROOT / "LICENSE", "LICENSE"),
                      (ROOT / "src/star_input/LICENSE", "STAR-input-LICENSE"),
                      (ROOT / "third_party/rapidmacs/LICENSE", "RapidMACS-LICENSE"),
                      (ROOT / "third_party/htslib/LICENSE", "HTSlib-LICENSE")):
        shutil.copy2(src, dest / "share/licenses/chromap-suite" / name)
    for src in (ROOT / "README.md", ROOT / "CHANGELOG.md", notes):
        shutil.copy2(src, dest / "share/doc/chromap-suite" / src.name)
    meta = dict(provenance(), release=tag, suite_version=suite, debian_upstream_version=deb_version,
                architecture=arch, distribution=distro,
                glibc=output(["getconf", "GNU_LIBC_VERSION"]).split()[-1],
                runtime_dependencies=runtime_dependencies(dest),
                binaries={b: sha256(dest / "bin" / b) for b in BINARIES})
    (dest / "share/chromap-suite/release.json").write_text(json.dumps(meta, indent=2) + "\n")
    (dest / "VERSION").write_text(f"chromap-suite {tag}\nsource-revision {meta['source_revision']}\n"
                                  f"platform linux-{arch}\nglibc-baseline {meta['glibc']}\n")
    return meta


def archive_tree(source, dest, prefix, epoch, exclude_debian=False):
    # Normalize timestamps and ownership; preserve executable modes.
    with dest.open("wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=epoch) as gz:
        with tarfile.open(fileobj=gz, mode="w") as archive:
            for path in [source, *sorted(source.rglob("*"))]:
                relative = path.relative_to(source)
                if exclude_debian and relative.parts and relative.parts[0] == "debian":
                    continue
                info = archive.gettarinfo(str(path), str(Path(prefix) / relative))
                info.uid = info.gid = 0
                info.uname = info.gname = ""
                info.mtime = epoch
                if info.isfile():
                    with path.open("rb") as data:
                        archive.addfile(info, data)
                else:
                    archive.addfile(info)


def build_tarball(stage, out):
    meta = json.loads((stage / "share/chromap-suite/release.json").read_text())
    name = f"Chromap-suite-{meta['release']}-linux-{meta['architecture']}-glibc{meta['glibc'].replace('.', '')}"
    path = out / f"{name}.tar.gz"
    archive_tree(stage, path, name, meta["source_date_epoch"])
    return path


def build_deb(stage, out):
    meta = json.loads((stage / "share/chromap-suite/release.json").read_text())
    version = f"{meta['debian_upstream_version']}-1~{meta['distribution']}.1"
    with tempfile.TemporaryDirectory() as temp:
        pkg = Path(temp) / "package"
        shutil.copytree(stage, pkg / "usr/lib/chromap-suite")
        (pkg / "usr/bin").mkdir(parents=True)
        for name in BINARIES:
            (pkg / "usr/bin" / name).symlink_to(f"../lib/chromap-suite/bin/{name}")
        doc = pkg / "usr/share/doc/chromap-suite"
        shutil.copytree(stage / "share/doc/chromap-suite", doc)
        shutil.copy2(ROOT / "debian/copyright", doc / "copyright")
        (pkg / "DEBIAN").mkdir()
        run(["dpkg-gencontrol", "-pchromap-suite", f"-v{version}", f"-P{pkg}",
             f"-O{pkg / 'DEBIAN/control'}", f"-f{Path(temp) / 'files'}",
             f"-Vshlibs:Depends={meta['runtime_dependencies']}", "-Vmisc:Depends="], cwd=ROOT)
        # GitHub release assets normalize '~'; keep filenames download-stable.
        path = out / f"chromap-suite_{version.replace('~', '.')}_{meta['architecture']}.deb"
        run(["dpkg-deb", "--root-owner-group", "--build", pkg, path],
            env=dict(os.environ, SOURCE_DATE_EPOCH=str(meta["source_date_epoch"])))
    return path


def export_source(dest, version):
    meta = dict(provenance(), release=version_info(version)[0])
    dest.mkdir(parents=True)
    for repo, revision, target in ((ROOT, meta["source_revision"], dest),
                                  (ROOT / "third_party/rapidmacs", meta["rapidmacs_revision"],
                                   dest / "third_party/rapidmacs")):
        target.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryFile() as archive:
            run(["git", "archive", "--format=tar", revision], cwd=repo, stdout=archive)
            archive.seek(0)
            run(["tar", "-xf", "-", "-C", target], stdin=archive)
    (dest / ".release-source.json").write_text(json.dumps(meta, indent=2) + "\n")
    return meta


def build_source(out, version):
    from email.utils import formatdate
    tag, suite, upstream = version_info(version)
    with tempfile.TemporaryDirectory() as temp:
        temp = Path(temp)
        source = temp / f"chromap-suite-{upstream}"
        meta = export_source(source, tag)
        archive_tree(source, temp / f"chromap-suite_{upstream}.orig.tar.gz", source.name,
                     meta["source_date_epoch"], exclude_debian=True)
        (source / "debian/changelog").write_text(
            f"chromap-suite ({upstream}-1) unstable; urgency=medium\n\n"
            f"  * Package Chromap Suite {tag}.\n\n"
            f" -- Ling-Hong Hung <lhhunghimself@gmail.com>  {formatdate(meta['source_date_epoch'], usegmt=True)}\n")
        run(["dpkg-source", "-b", source], cwd=temp,
            env=dict(os.environ, SOURCE_DATE_EPOCH=str(meta["source_date_epoch"])))
        paths = [p for p in temp.iterdir() if p.is_file()]
        for path in paths:
            shutil.copy2(path, out / path.name)
        bundle = temp / "debian-source"
        bundle.mkdir()
        for path in paths:
            shutil.copy2(path, bundle / path.name)
        archive_tree(bundle, out / f"Chromap-suite-{tag}-debian-source.tar.gz",
                     bundle.name, meta["source_date_epoch"])
    return [out / p.name for p in paths]


def build(args):
    tag, suite, _ = version_info(args.version)
    out = args.out_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if any(out.iterdir()):
        raise ValueError(f"use an empty output directory to avoid stale artifacts: {out}")
    tests = args.test_dir.resolve()
    tests.mkdir(parents=True, exist_ok=True)
    # The packaging operations below cannot run if compilation or any test fails.
    with (tests / "build.log").open("w") as log:
        run(["make", f"-j{args.jobs}", "all", "chromap_lib_runner", "tests/fastq_intake_harness"],
            cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    run(["bash", ROOT / "scripts/release/run_release_tests.sh", tests], cwd=ROOT)
    with tempfile.TemporaryDirectory(dir=out.parent) as temp:
        stage = Path(temp) / "stage"
        meta = stage_release(stage, tag)
        run(["bash", ROOT / "tests/run_release_artifact_smoke.sh", stage, suite,
             meta["source_revision"], tests / "staged-runtime"])
        run(["bash", ROOT / "tests/run_release_sdk_smoke.sh", stage, tests / "sdk"])
        tarball = build_tarball(stage, out)
        deb = build_deb(stage, out)
        (out / "build-manifest.json").write_text(json.dumps(dict(meta, tests="PASS",
            artifacts={p.name: sha256(p) for p in (tarball, deb)}), indent=2) + "\n")
    run(["bash", ROOT / "scripts/release/create_checksums.sh", out])
    print(f"PASS: built tested release artifacts in {out}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    build_parser = sub.add_parser("build", help="build, test, stage, check SDK and package tarball + deb")
    build_parser.add_argument("--version", required=True)
    build_parser.add_argument("--out-dir", type=Path, required=True)
    build_parser.add_argument("--test-dir", type=Path, default=ROOT / "plans/artifacts/release")
    build_parser.add_argument("--jobs", type=int, default=4)
    for command in ("stage", "source", "snapshot", "tarball"):
        p = sub.add_parser(command)
        p.add_argument("--version", default=provenance().get("release", "v" + suite_version()))
        p.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    os.chdir(ROOT)
    if args.command == "build":
        if args.jobs < 1:
            parser.error("--jobs must be positive")
        build(args)
    elif args.command == "stage":
        stage_release(args.out_dir.resolve(), args.version)
    elif args.command == "snapshot":
        export_source(args.out_dir.resolve(), args.version)
    else:
        args.out_dir.mkdir(parents=True, exist_ok=True)
        if args.command == "source":
            build_source(args.out_dir.resolve(), args.version)
        else:
            with tempfile.TemporaryDirectory() as temp:
                stage = Path(temp) / "stage"
                stage_release(stage, args.version)
                print(build_tarball(stage, args.out_dir.resolve()))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, subprocess.CalledProcessError) as error:
        raise SystemExit(f"release: {error}")
