#!/usr/bin/env python3
"""Generate a Fortran dispatcher for all loop_N subroutines in thermal_lines.f90."""

from __future__ import annotations

import re
from pathlib import Path

THERMAL_LINES = Path("thermal_lines.f90")
OUTPUT = Path("operator_dispatch.f90")

SUBROUTINE_RE = re.compile(
    r"^\s*subroutine\s+loop_(\d+)\s*\(([^)]*)\)",
    re.IGNORECASE | re.MULTILINE,
)


def normalise_arguments(argument_text: str) -> list[str]:
    return [argument.strip().lower() for argument in argument_text.split(",") if argument.strip()]


def main() -> None:
    source = THERMAL_LINES.read_text()
    loops: dict[int, list[str]] = {}

    for match in SUBROUTINE_RE.finditer(source):
        loop_number = int(match.group(1))
        arguments = normalise_arguments(match.group(2))
        loops[loop_number] = arguments

    if not loops:
        raise SystemExit("No loop_N subroutines found in thermal_lines.f90")

    loop_numbers = sorted(loops)

    lines: list[str] = []
    lines.append("module operator_dispatch")
    lines.append("    use iso_fortran_env, only : real64")
    lines.append("    use parameters")
    lines.append("    use lattice")
    lines.append("    use thermal_lines")
    lines.append("    implicit none")
    lines.append("")
    lines.append(f"    integer, parameter :: NUMBER_OF_LOOPS = {len(loop_numbers)}")
    lines.append("    integer, parameter :: LOOP_NUMBERS(NUMBER_OF_LOOPS) = [ &")

    for index, loop_number in enumerate(loop_numbers):
        separator = ", &" if index < len(loop_numbers) - 1 else " &"
        lines.append(f"        {loop_number}{separator}")

    lines.append("    ]")
    lines.append("")
    lines.append("contains")
    lines.append("")
    lines.append("    subroutine call_loop(loop_number, in_site, out_site, A11, blocking_level, gauge_field_blocked)")
    lines.append("        integer, intent(in) :: loop_number, in_site, blocking_level")
    lines.append("        integer, intent(out) :: out_site")
    lines.append("        complex(real64), intent(inout) :: A11(NCOL, NCOL)")
    lines.append("        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)")
    lines.append("")
    lines.append("        select case(loop_number)")

    for loop_number in loop_numbers:
        arguments = loops[loop_number]
        lines.append(f"        case({loop_number})")
        if "gauge_field_blocked" in arguments:
            lines.append(
                f"            call loop_{loop_number}(in_site, out_site, A11, blocking_level, gauge_field_blocked)"
            )
        else:
            lines.append(
                f"            call loop_{loop_number}(in_site, out_site, A11, blocking_level)"
            )

    lines.append("        case default")
    lines.append("            error stop 'Invalid loop number in generated call_loop'")
    lines.append("        end select")
    lines.append("    end subroutine call_loop")
    lines.append("")
    lines.append("end module operator_dispatch")
    lines.append("")

    OUTPUT.write_text("\n".join(lines))
    print(f"Generated {OUTPUT} with {len(loop_numbers)} loop dispatch cases")


if __name__ == "__main__":
    main()
