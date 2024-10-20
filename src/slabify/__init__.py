"""CLI init function."""
# These imports are necessary to register CLI commands. Do not remove!
from .boundary_mask_auto import create_boundary_mask_auto  # noqa: F401
from .slabify import slabify  # noqa: F401
from .stopgap_tm import sg_tm_create_boundary_mask  # noqa: F401
