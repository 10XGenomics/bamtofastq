workspace(name = "bamtofastq")

load(
    "@bazel_tools//tools/build_defs/repo:git.bzl",
    "git_repository",
)

# When debugging iteratively, it's useful to switch to a local repo rather
# than repeated pushes to github.  To do this, run with
# --override_repository=tenx_bazel_rules=/path/to/repository
git_repository(
    name = "tenx_bazel_rules",
    commit = "69af1bc334dc631d16c928af61ef0e27d8157d91",
    remote = "https://github.com/10XDev/tenx_bazel_rules.git",
)

load(
    "@tenx_bazel_rules//:deps.bzl",
    "toolchain_dependencies",
)

toolchain_dependencies()

load(
    "@tenx_bazel_rules//:toolchains.bzl",
    "register_pipeline_toolchains",
)

register_pipeline_toolchains()

load(
    "@tenx_bazel_rules//rules:cargo_repository.bzl",
    "cargo_repository",
)

cargo_repository(
    name = "bamtofastq_cargo_dependencies",
    lockfile = ":Cargo.lock",
    srcs = [
        ":Cargo.toml",
    ],
)
