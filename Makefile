
.PHONY: all fmt fmt-check lint test test-ignored build release clean check

all: fmt lint test build

fmt:
	cargo fmt --all

fmt-check:
	cargo fmt --all -- --check

lint:
	cargo clippy --workspace -- -D warnings

test:
	cargo test --workspace

# Run #[ignore]-gated concordance harnesses (GRCh38 + GRCh37). Requires both
# assemblies' data files on disk (see `vareffect setup --assembly all`).
test-ignored:
	FASTA_PATH=data/vareffect/GRCh38.bin \
	GRCH37_FASTA=data/vareffect/GRCh37.bin \
	GRCH37_TRANSCRIPTS=data/vareffect/transcript_models_grch37.bin \
	  cargo test -p vareffect --release -- --ignored --nocapture

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