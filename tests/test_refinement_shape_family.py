"""``flexref``, ``emref`` and ``mdref`` are one job shape.

They differ in their CNS recipe and in nothing else that a schedule can see:
each takes one input model, refines it, and writes one model out.  Test sets
that cover job shapes -- including the caching contract suite -- therefore
cover the shape through a single representative, and adding the same case
again for the other two would buy nothing.

That is sound coverage only while the premise holds, so the premise is a
constraint on the code:

    a change to the job emission, seeding, indexing or cache participation of
    one of the three must be applied to all three, or the divergence must be
    deliberate and written down.

Violating it does not turn any test red.  The suite stays green while two of
the three modules carry the defect, because the representative has quietly
stopped representing -- which is worse than a coverage gap, since a gap is at
least visible in the case list.

These tests pin the part of that constraint which can be checked cheaply: the
three modules schedule their replicas and derive their seeds through one
shared implementation rather than through three copies of it.
"""

import pytest

from haddock.libs import libcns
from haddock.modules.refinement import emref, flexref, mdref


REFINEMENT_FAMILY = (flexref, emref, mdref)


@pytest.mark.parametrize("module", REFINEMENT_FAMILY, ids=lambda m: m.__name__)
def test_refinement_family_shares_one_replica_schedule(module):
    """One emission order for the shape, not one per module."""
    assert module.refinement_schedule is libcns.refinement_schedule


@pytest.mark.parametrize("module", REFINEMENT_FAMILY, ids=lambda m: m.__name__)
def test_refinement_family_shares_one_seed_derivation(module):
    """One seeding rule for the shape, not one per module."""
    assert module.derive_seed is libcns.derive_seed
