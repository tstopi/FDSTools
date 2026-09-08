# -*- coding: utf-8 -*-
"""Pipeline exception types (kept separate to avoid circular imports)."""


class PipelineError(RuntimeError):
    """Base class for any pipeline phase failure."""


class CollectError(PipelineError):
    """Raised when source geometry cannot be gathered (Phase 2)."""
