# Doping-Calculator
A calculator for calculating the mass of non-stoichiometric compound
Make sure that .py file and relative atomic mass file are in the same folder
用法示例/Usage:
```
(base) storm@DESKTOP-HE4FQ8Q:~/my-learn/trash-vaspkit/test/test$ ls
Doping-Calculator.py  atomic_masses.txt
(base) storm@DESKTOP-HE4FQ8Q:~/my-learn/trash-vaspkit/test/test$ python Doping-Calculator.py
High-Precision Compound Synthesis Calculator
============================================================
[Success] Atomic mass data loaded. Elements: Li, P, S, Cl, Br, Na, O, F, I, Sc, Ti
Enter target compound formula (e.g. Li5.3PS4.3Cl0.7Br): Li5.5PS4.5Cl1.5
Parsed target compound: {'Li': Decimal('5.5'), 'P': Decimal('1'), 'S': Decimal('4.5'), 'Cl': Decimal('1.5')}
[Calculated] Target molar mass = 266.599 g/mol
Enter number of reagents: 3
Reagent 1 formula: LiCl
[Calculated] LiCl molar mass = 42.394 g/mol
Reagent 2 formula: P2S5
[Calculated] P2S5 molar mass = 222.248 g/mol
Reagent 3 formula: Li2S
[Calculated] Li2S molar mass = 45.942 g/mol
Enter input type for target (0 for mass in g, 1 for moles in mol): 0
Enter the target amount: 1
[Calculated] Target amount in moles = 0.00375095180402027014354892553985573839361738041027910832373714830138147554942066549386906927632886845 mol

============================================================
To synthesize 1 g of target compound:
Reagent 1 (LiCl):
  Moles required: 0.0056264277 mol
  Mass required:  0.2385267762 g

Reagent 2 (P2S5):
  Moles required: 0.0018754759 mol
  Mass required:  0.4168207683 g

Reagent 3 (Li2S):
  Moles required: 0.0075019036 mol
  Mass required:  0.3446524556 g

Total reagent mass: 1.0000000000 g

Verification of element amounts:
Cl: Target = 0.0056264277, Actual = 0.0056264277, Diff = -5.215323e-18
Li: Target = 0.0206302349, Actual = 0.0206302349, Diff = -5.789519e-18
P: Target = 0.0037509518, Actual = 0.0037509518, Diff = -1.014355e-17
S: Target = 0.0168792831, Actual = 0.0168792831, Diff = -2.564597e-17
```
