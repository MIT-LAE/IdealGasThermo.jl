# Deprecation machinery for the legacy mutable layer (`Gas{N}`, `Gas1D` and the
# Dict-combustion / mutable-turbo functions they back). See ADR-0002 (the layer is
# superseded by FrozenGas/GasState) and ADR-0007 (the v2.0 migration ships these
# loudly deprecated in a `2.0.0-betaN` series, then deletes them in 2.0.0 final).
#
# Policy (ADR-0007): the warning is LOUD (a plain `@warn`, always shown — not
# `Base.depwarn`, which Julia hides unless `--depwarn=yes`) but fires at most once
# per entry point per session (`maxlog = 1`, keyed by `_id`). The only runtime
# choke points are the `Gas` and `Gas1D` *constructors*: every legacy public
# function (`set_TP!`, `print_thermo_table`, `gas_burn`, …) takes a `Gas`/`Gas1D`,
# so constructing one is the single signal that a caller is on the legacy path.
# Those downstream functions carry docstring deprecation notes but no runtime
# warning, which would otherwise cascade through internal plumbing.

"""
    _legacy_warn(name::Symbol, replacement::AbstractString)

Emit a one-time (`maxlog = 1`, per `name`) deprecation warning for a legacy
entry point. `name` is the deprecated symbol (e.g. `:Gas`); `replacement` names
the pure-core API to migrate to (e.g. `"FrozenGas / GasState"`).
"""
function _legacy_warn(name::Symbol, replacement::AbstractString)
    @warn "`$(name)` is deprecated and will be removed in IdealGasThermo v2.0.0; " *
          "use $(replacement) instead. (ADR-0002, ADR-0007)" maxlog = 1 _id = name
    return nothing
end
