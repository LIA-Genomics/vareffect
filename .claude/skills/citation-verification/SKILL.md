---
name: citation-verification
description: >-
  Anti-fabrication discipline: verify every PMID/DOI, identifier, threshold, version, API/type/field name or guideline statement against a fetched source or opened file, or drop it. Use whenever a task asserts or cites something that must be true.
---

# Citation & claim verification — never assert from memory

Your credibility rests on a single rule: **every factual claim traces to something you
actually verified, or it does not ship.** A confident, wrong claim is worse than no claim
— it corrupts everything downstream. Precision over recall.

## Prime directive — zero fabrication

Every PMID/DOI, author, title, year, journal, gene/identifier, threshold, version number,
and type/field/function/API name you output MUST trace to **a source you actually fetched**
or **a location you actually opened** (cite it — a URL, or a `path:line`). Never model
memory. If you cannot ground a claim, you **drop it** and say so. When you are unsure
whether a name/number is real, treat it as unverified, not as true.

## The verify-loop for an identifier (PMID/DOI/version/etc.)

1. **Format ≠ existence.** A string that is *well-formed* (e.g. a 1–8-digit PMID, a valid
   semver) is not thereby *real*. Well-formedness is a cheap pre-check, never proof.
2. **Fetch it.** Resolve the identifier against its authoritative registry (e.g. NCBI
   E-utilities `efetch` for a PMID, the package registry for a crate version, the real
   source file for a symbol). If it does not resolve → **drop it**.
3. **Content must match the claim.** A **real identifier attached to the wrong thing is
   still fabrication.** A PMID that exists but points at an unrelated paper does not support
   your statement — verify that the fetched source actually says what you claim. (Illustrative
   failures: a real PMID for an LDL GWAS, and another for an orthopedic note, each wrongly
   cited into a disease-prevalence claim — both real IDs, both fabrications in context.)

## The truncation / empty-input guard

If the input you were given looks **empty, truncated, or undefined**, STOP and flag it —
do **not** backfill from memory to "complete" it. Inventing plausible-looking content to
fill a gap is the most dangerous failure mode (e.g. hallucinating a list of items because
the real input didn't arrive). A flagged gap is a correct result; a confident fabrication
is not.

## Authoritative-source hierarchy

Prefer **primary, authoritative** sources; never a blog, forum, marketing page, or preprint
as the *sole* anchor for a consequential claim.
- **Variant annotation (vareffect's domain):** the HGVS Nomenclature recommendations
  (hgvs-nomenclature.org), the Sequence Ontology, Ensembl VEP documentation and the
  `ensembl-vep` / `ensembl-variation` source at the release in question, the Ensembl REST
  API's actual response, NCBI (RefSeq, MANE, genetic code tables, the GFF3 and GTF
  documentation for RefSeq), GENCODE's format documentation, the GFF3 specification
  (Sequence Ontology), the VCF specification (samtools/hts-specs), and for downstream
  clinical use ClinGen SVI publications on PubMed. A secondary source may corroborate, but
  the anchor is primary.
- **Software:** the real source file / the official API docs / the package registry — not a
  Q&A site's recollection.
- **Behavior beats prose.** When VEP's documentation and its observed output disagree, record
  both; the observed output (request and response) is what concordance is measured against.

## Current behavior and required behavior

Code, tests and observations establish current implementation behavior. The binding charter
and authorized requirements establish required behavior. A disagreement can be an implementation
defect; do not automatically classify it as documentation drift. Descriptive docs and memory are
orientation to verify against source, while normative obligations retain their authority.

## Output discipline

- **Precision over recall — never pad.** An ungrounded finding is noise. Ship fewer,
  verified claims rather than more, hedged ones.
- Separate what you **verified** from what you **could not** — put unverifiable items in an
  explicit "open questions / unverifiable" bucket, never silently among your assertions.
- Cite the anchor inline for each non-obvious claim (URL or `path:line`), so a reader can
  re-check without trusting you.
