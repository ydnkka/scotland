"""Compatibility imports for the renamed sensitivity figure module."""

from .sensitivity_figures import (
    make_clade_association_figures,
    make_clade_association_outputs,
    make_clade_association_summary_tables,
)

__all__ = [
    "make_clade_association_figures",
    "make_clade_association_outputs",
    "make_clade_association_summary_tables",
]
