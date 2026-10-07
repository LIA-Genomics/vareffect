
.PHONY: all fmt fmt-check lint test doc-check build release clean check install-cli harness-check

all: harness-check fmt lint test doc-check build

fmt:
	cargo fmt --all

fmt-check:
	cargo fmt --all -- --check

lint:
	cargo clippy --workspace --all-targets -- -D warnings

test:
	cargo test --workspace

doc-check:
	RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps

build:
	cargo build --workspace

release:
	cargo build --release --workspace

clean:
	cargo clean

check:
	cargo check --workspace

install-cli:
	cargo install --path vareffect-cli --force

harness-check:
	python3 scripts/vareffect_harness.py check
	python3 -B scripts/test_vareffect_harness.py
