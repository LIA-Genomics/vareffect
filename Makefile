
.PHONY: all fmt fmt-check lint test build release clean check

all: fmt lint test build

fmt:
	cargo fmt --all

fmt-check:
	cargo fmt --all -- --check

lint:
	cargo clippy --workspace -- -D warnings

test:
	cargo test --workspace

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