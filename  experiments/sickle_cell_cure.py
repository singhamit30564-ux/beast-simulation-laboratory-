"""
Sickle Cell Disease Cure Simulation
Target: BCL11A erythroid enhancer disruption
Strategy: Reactivate fetal hemoglobin (HbF)
"""

from src.crispr.cas9_sim import Cas9Simulator
from src.genome.pam_scanner import PAMScanner

def simulate_sickle_cell_cure():
    """Simulate CRISPR therapy for sickle cell disease."""
    
    # BCL11A enhancer region (simplified sequence)
    bcl11a_enhancer = (
        "ATCGATCGATCGATCGATCGATCGATCGATCG"
        "GAGTCCATGGTGGCCATGTGGTGGCCATGTGG"
        "TCCATGGTGGCCATGTGGTGGCCATGTGGTGG"
        "CCATGTGGTGGCCATGTGGTGGCCATGTGGCC"
    )
    
    print("=" * 50)
    print("SICKLE CELL DISEASE CURE SIMULATION")
    print("Target: BCL11A erythroid enhancer")
    print("Strategy: Disrupt enhancer → Reactivate HbF")
    print("=" * 50)
    
    # Initialize Cas9
    cas9 = Cas9Simulator(variant="SpCas9")
    
    # Find PAM sites
    pams = PAMScanner.find(bcl11a_enhancer, variant="SpCas9")
    print(f"\nFound {len(pams)} PAM sites in target region")
    
    if not pams:
        print("No valid PAM sites found!")
        return
    
    # Select best PAM (first one for demo)
    target_pam = pams[0]
    print(f"\nSelected PAM: {target_pam.sequence} at position {target_pam.position}")
    
    # Simulate cut with NHEJ (creates indel = enhancer disruption)
    result = cas9.cut(
        target=bcl11a_enhancer,
        pam=target_pam,
        repair_pathway="NHEJ",
        cell_type="iPSC"  # Patient-derived stem cells
    )
    
    print(f"\n--- RESULTS ---")
    print(f"Cut Efficiency: {result.efficiency:.1f}%")
    print(f"Repair Pathway: {result.repair_pathway}")
    print(f"Cut Position: {result.cut_position}")
    
    print(f"\nPredicted Indel Profile:")
    for indel, prob in result.indels.items():
        print(f"  {indel}: {prob}%")
    
    print(f"\nOff-targets detected: {len(result.off_targets)}")
    for ot in result.off_targets:
        print(f"  - {ot['chromosome']}:{ot['position']} "
              f"(mismatches: {ot['mismatches']}, CFD: {ot['cfd_score']})")
    
    print(f"\nClinical Outcome Prediction:")
    if result.efficiency > 70 and len(result.off_targets) < 3:
        print("  ✅ High efficacy predicted. HbF reactivation likely.")
        print("  ✅ Potential functional cure for sickle cell disease.")
    else:
        print("  ⚠️  Optimization needed. Consider alternative gRNA.")

if __name__ == "__main__":
    simulate_sickle_cell_cure()
