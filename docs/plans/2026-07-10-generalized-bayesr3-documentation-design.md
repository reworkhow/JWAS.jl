# Generalized BayesR3 Documentation Design

## Goal

Reframe the JWAS fast-block documentation around a generalized BayesR3 strategy.
The documentation should credit the BayesR3 blocked Gibbs design while making
clear that JWAS deliberately supports additional statistical models, user-defined
unequal blocks, and multiple execution schedules.

This is a documentation-only change. It does not change sampler behavior or the
current `fast_blocks` API.

## Terminology

The **BayesR3 block strategy** consists of three computational steps:

1. construct the right-hand side for the current marker block;
2. perform repeated conditional marker updates using the block Gram matrix;
3. reconcile the block effect change with the global residual when leaving the
   block.

JWAS implements a **generalized BayesR3 strategy** because it reuses this
computational structure beyond the original BayesR mixture model and allows the
block partition and execution schedule to vary.

The documentation must not claim that JWAS reproduces the original BayesR3
software or every detail of the paper's scheduling convention.

## Independent Design Dimensions

The documentation will separate three choices that are currently easy to
conflate:

1. **Statistical model**
   - BayesA/B/C
   - BayesR
   - supported annotated and multi-trait variants
2. **Block partition**
   - automatically chosen block size
   - user-supplied fixed block size
   - user-supplied ordered block starts, which may produce unequal contiguous
     blocks
3. **Execution schedule**
   - exact sequential blocks using the current residual at each block boundary
   - approximate independent blocks using a frozen sweep-level residual, except
     in the special case of vanishing off-block weighted crossproducts

## User-Defined Blocks

JWAS users may define block boundaries from scientific knowledge such as LD,
recombination, IBD, pedigree structure, or genomic regions. Blocks are not
required to be equal in size.

The current API accepts ordered marker start positions, so the resulting blocks
are contiguous in the current marker order. Arbitrary non-contiguous marker
membership is outside this documentation change and may be considered as a
future API feature.

The documentation may state that an informed partition can improve computational
efficiency or mixing, especially when strongly coupled markers are grouped
together. It must not guarantee better convergence: partition quality, block
size, repeat policy, and the posterior correlation structure all affect mixing.

## Block and Repeat Notation

The revised page will distinguish:

- `b`: nominal block size;
- `s_i`: realized size of block `i`;
- `r_i`: number of inner repeats applied to block `i`.

For the original BayesR3 convention, `r_i = b` for every block, including a
final block with `s_i < b`. Under the paper's chain-length convention, this gives
every SNP effect the same number of updates.

In the current JWAS BayesC block path, `r_i = s_i`. The current BayesR block path
uses one repeat during burn-in and `r_i = s_i` afterward. This distinction must
be stated without implying that different repeat policies change the underlying
block residual identity.

## Exactness and Approximation

Sequential block mode is a different transition schedule from a dense
marker-by-marker sweep, but it composes conditional Gibbs updates and refreshes
the global residual before processing the next block. It should be described as
the exact sequential block schedule for the implemented model, not as a
heuristic approximation to dense Gibbs.

`independent_blocks=true` updates blocks from a frozen sweep-level residual. It
is the explicitly approximate mode unless off-block weighted genotype
crossproducts vanish.

## Documentation Changes

The main page will remain at `docs/src/manual/block_bayesc.md` to preserve
existing links, but its heading and opening will describe the BayesR3 block
strategy in JWAS.

The revision will:

- replace the broad statement that JWAS uses BayesC rather than BayesR;
- describe the worked BayesC comparison as one application of the generalized
  strategy;
- add a clear section on user-defined unequal blocks;
- separate model, partition, repeat policy, and execution schedule;
- qualify convergence claims;
- distinguish the JWAS operation-count model from the BayesR3 empirical runtime
  fit and avoid presenting them as a software benchmark;
- retain links to the original BayesR3 paper;
- align related wording in the Annotated BayesR, workflow, and navigation text
  where needed.

## Verification

After editing:

1. search the documentation for stale claims that sequential block mode or
   BayesR3 is merely an approximation;
2. inspect the rendered wording for consistent definitions of `b`, `s_i`, and
   `r_i`;
3. run `julia --project=docs --startup-file=no docs/make.jl`;
4. record the completed changes and verification result in a matching
   implementation note under `docs/plans/`.
