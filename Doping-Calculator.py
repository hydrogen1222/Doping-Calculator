from decimal import Decimal, getcontext
import re
import os
from sympy import symbols, Eq, nsimplify, Matrix, linsolve
from collections import defaultdict

# Set high precision for Decimal calculations
getcontext().prec = 100

class FormulaParser:
    @staticmethod
    def parse(formula: str) -> dict:
        """
        Parse a chemical formula and return a dictionary of element counts.
        For example, "Li5.5PS4.5Cl1.5" → {'Li': Decimal('5.5'), 'P': Decimal('1'), 'S': Decimal('4.5'), 'Cl': Decimal('1.5')}
        """
        elements = {}
        # Regex pattern: captures an uppercase letter followed by optional lowercase letter (element symbol),
        # then an optional number (integer or decimal) for its coefficient.
        pattern = re.compile(r"([A-Z][a-z]?)([\d.]*)?")
        index = 0

        while index < len(formula):
            match = pattern.match(formula, index)
            if not match:
                raise ValueError(f"Invalid formula format: {formula[index:]}")

            elem = match.group(1)               # Element symbol (e.g., "Li", "P", "Cl")
            number = match.group(2)             # Coefficient string (e.g., "5.5" or "")

            # Default coefficient = 1 if no number is provided
            if number == '' or number is None:
                coeff = Decimal(1)
            else:
                try:
                    coeff = Decimal(number)
                except:
                    raise ValueError(f"Invalid coefficient: {number}")

            # Accumulate coefficient if element appears multiple times
            if elem in elements:
                elements[elem] += coeff
            else:
                elements[elem] = coeff

            index = match.end()

        return elements

class AtomicMassReader:
    @staticmethod
    def read(filename="atomic_masses.txt") -> dict:
        """
        Read an atomic mass file with lines like "Li 6.941" and return a dict {element: Decimal(mass)}.
        Raises FileNotFoundError if the file does not exist or ValueError if contents are malformed.
        """
        atomic_masses = {}
        line_num = 0

        try:
            with open(filename, 'r') as f:
                for line in f:
                    line_num += 1
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue  # Skip blank lines and comments

                    parts = line.split()
                    if len(parts) != 2:
                        raise ValueError(f"Line {line_num} format error, expected 'element mass'")

                    element, mass = parts[0], parts[1]
                    try:
                        atomic_masses[element] = Decimal(mass)
                    except:
                        raise ValueError(f"Line {line_num} invalid mass value: {mass}")

        except FileNotFoundError:
            raise FileNotFoundError(f"Atomic mass file {filename} not found")

        if not atomic_masses:
            raise ValueError("Atomic mass file is empty")

        return atomic_masses

class SynthesisCalculator:
    def __init__(self, atomic_masses: dict):
        """
        Initialize the calculator with a dictionary of atomic masses (element → Decimal mass).
        """
        self.target = {}               # Dict of element counts for target compound
        self.target_mm = Decimal(0)    # Molar mass of the target compound
        self.reagents = []             # List of dicts for each reagent's element counts
        self.reagent_formulas = []     # Original string formulas for reagents
        self.molar_masses = []         # Computed molar mass (Decimal) for each reagent
        self.atomic_masses = atomic_masses

    def input_target(self) -> Decimal:
        """
        Prompt user to enter the target compound formula, parse it, calculate and print its molar mass.
        Returns the computed molar mass (Decimal).
        """
        formula = input("Enter target compound formula (e.g. Li5.3PS4.3Cl0.7Br): ").strip()
        self.target = FormulaParser.parse(formula)

        # Calculate target molar mass by summing (element count × atomic mass)
        for elem, coeff in self.target.items():
            if elem not in self.atomic_masses:
                raise ValueError(f"Element {elem} mass not defined")
            self.target_mm += coeff * self.atomic_masses[elem]

        print(f"Parsed target compound: {self.target}")
        print(f"[Calculated] Target molar mass = {self.target_mm.normalize()} g/mol")
        return self.target_mm

    def input_reagents(self):
        """
        Prompt user for number of reagents, then for each reagent's formula.
        Parse each reagent, compute its molar mass, and store results.
        """
        num = int(input("Enter number of reagents: "))
        for i in range(num):
            formula = input(f"Reagent {i+1} formula: ").strip()
            self.reagent_formulas.append(formula)
            try:
                parsed = FormulaParser.parse(formula)
                mm = Decimal(0)
                for elem, coeff in parsed.items():
                    if elem not in self.atomic_masses:
                        raise ValueError(f"Element {elem} mass not defined")
                    mm += coeff * self.atomic_masses[elem]

                self.reagents.append(parsed)
                self.molar_masses.append(mm)
                print(f"[Calculated] {formula} molar mass = {mm.normalize()} g/mol")

            except Exception as e:
                print(f"Input error: {e}")
                exit()

    def validate_elements(self):
        """
        Ensure that all elements in the target appear in at least one reagent.
        Warn if reagents contain elements not in the target.
        """
        target_elements = set(self.target.keys())
        provided = set()
        for r in self.reagents:
            provided.update(r.keys())

        missing = target_elements - provided
        extra = provided - target_elements

        if missing:
            print(f"Error: Missing elements: {', '.join(missing)}")
            exit()
        if extra:
            print(f"Warning: Extra elements provided: {', '.join(extra)}")

    def build_equations(self, n: Decimal):
        """
        Build the coefficient matrix A and vector b for the linear system A x = b,
        where x = [x0, x1, ..., x_{num_reagents - 1}] are moles of each reagent,
        and n is the number of moles of target compound.
        Returns:
            A: list of lists (coefficients of each element in each reagent)
            b: list (required total moles of each element = target_coeff × n)
            all_elements: sorted list of all elements involved
        """
        # Collect all elements present in either target or any reagent
        all_elements = set(self.target.keys())
        for reagent in self.reagents:
            all_elements.update(reagent.keys())
        all_elements = sorted(all_elements)

        A = []
        b = []
        for element in all_elements:
            row = []
            for reagent in self.reagents:
                coeff = reagent.get(element, Decimal(0))
                row.append(coeff)
            A.append(row)

            target_coeff = self.target.get(element, Decimal(0))
            b.append(target_coeff * n)

        return A, b, all_elements

    def solve(self, n: Decimal):
        """
        Solve the linear system for reagent mole amounts x0, x1, ..., x_{m-1},
        where m = number of reagents, subject to A x = b.
        If no exact solution exists, print a clear message and exit.
        Returns:
            solution_dicts: a list of dictionaries mapping each Symbol xi to a Decimal value.
            elements: the list of all elements used in the system (for verification).
        """
        A, b, elements = self.build_equations(n)
        num_reagents = len(self.reagents)

        # Convert Decimal entries to float for Sympy solving
        A_float = [[float(val) for val in row] for row in A]
        b_float = [float(val) for val in b]

        # Convert lists to Sympy matrices for linsolve
        A_mat = Matrix(A_float)
        b_mat = Matrix(b_float)

        # Define symbolic variables x0, x1, ..., x_{num_reagents-1}
        vars = symbols(f'x0:{num_reagents}')

        # Attempt to solve A x = b exactly (with floats that represent exact .5 increments correctly)
        solution = linsolve((A_mat, b_mat), vars)

        # If linsolve returns an empty set (no exact solution), prompt and exit
        if len(solution) == 0:
            print("Error: No exact solution exists. Cannot synthesize the target composition with the given reagents.")
            exit()

        # Convert the exact solution(s) in Sympy FiniteSet to list of dicts
        sol_list = list(solution)  # Each element is a tuple of values for x0, x1, ...
        solution_dicts = []
        for sol_tuple in sol_list:
            sol_dict = {}
            for i, val in enumerate(sol_tuple):
                # Use nsimplify to clean up any small floating-point noise, then convert to string → Decimal
                cleaned_val = nsimplify(val).evalf(100)
                sol_dict[vars[i]] = Decimal(str(cleaned_val))
            solution_dicts.append(sol_dict)

        return solution_dicts, elements

    def calculate(self):
        """
        Main routine to prompt for target amount (mass or moles), solve for reagent amounts,
        and print out the required moles and masses for each reagent. Also verifies element balances.
        """
        input_type = input("Enter input type for target (0 for mass in g, 1 for moles in mol): ").strip()
        if input_type not in ['0', '1']:
            print("Error: Invalid input type. Please enter 0 or 1.")
            exit()

        value = input("Enter the target amount: ").strip()
        try:
            value = Decimal(value)
        except:
            print(f"Invalid amount: {value}")
            exit()

        # Convert mass → moles if needed
        if input_type == '0':
            n = value / self.target_mm
            print(f"[Calculated] Target amount in moles = {n.normalize()} mol")
        else:
            n = value

        # Solve for reagent moles (this will exit if no exact solution)
        solutions, elements = self.solve(n)
        solution = solutions[0]  # Take the first exact solution

        print("\n" + "="*60)
        print(f"To synthesize {value} {'g' if input_type == '0' else 'mol'} of target compound:")

        total_mass = Decimal(0)
        for i, var in enumerate(solution.keys()):
            moles = solution[var]
            if moles < 0:
                print(f"Error: Negative moles for reagent {i+1}")
                exit()

            mass = moles * self.molar_masses[i]
            total_mass += mass

            print(f"Reagent {i+1} ({self.reagent_formulas[i]}):")
            print(f"  Moles required: {moles.normalize():.10f} mol")
            print(f"  Mass required:  {mass.normalize():.10f} g\n")

        print(f"Total reagent mass: {total_mass.normalize():.10f} g")

        # Verification: check that elemental totals match target composition
        print("\nVerification of element amounts:")
        element_totals = defaultdict(Decimal)
        for i, reagent in enumerate(self.reagents):
            moles = solution[symbols(f'x{i}')]
            for element, coeff in reagent.items():
                element_totals[element] += moles * coeff

        for element in elements:
            target_amount = self.target.get(element, Decimal(0)) * n
            actual_amount = element_totals.get(element, Decimal(0))
            diff = (actual_amount - target_amount).normalize()
            print(f"{element}: Target = {target_amount:.10f}, Actual = {actual_amount:.10f}, Diff = {diff:.6e}")

def main():
    print("High-Precision Compound Synthesis Calculator")
    print("="*60)

    # Read atomic masses from file
    try:
        atomic_masses = AtomicMassReader.read()
        print(f"[Success] Atomic mass data loaded. Elements: {', '.join(atomic_masses.keys())}")
    except Exception as e:
        print(f"Atomic mass file error: {e}")
        print("Ensure atomic_masses.txt exists in current directory with lines like:")
        print("Li 6.941")
        print("P 30.973762")
        print("S 32.065")
        print("Cl 35.453")
        print("Br 79.904")
        exit()

    calc = SynthesisCalculator(atomic_masses)
    calc.input_target()
    calc.input_reagents()
    calc.validate_elements()
    calc.calculate()

if __name__ == "__main__":
    main()
