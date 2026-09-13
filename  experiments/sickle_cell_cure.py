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
    print("🩸 SICKLE CELL DISEASE CURE SIMULATION 🩸")
    print("Target: BCL11A erythroid enhancer")
    print("Strategy: Disrupt enhancer → Reactivate HbF")
    print("=" * 50)

    # Initialize Cas9 with a seed for reproducible results
    cas9 = Cas9Simulator(variant="SpCas9", seed=42)
    
    # Initialize Scanner
    scanner = PAMScanner()

    # Find PAM sites
    pams = scanner.find(bcl11a_enhancer, variant="SpCas9")
    print(f"\n🔍 Found {len(pams)} PAM sites in target region")

    if not pams:
        print("❌ No valid PAM sites found! Exiting.")
        return

    # Smart PAM Selection: Prefer a PAM that has at least 20bp upstream for a full guide RNA
    valid_pams = [p for p in pams if p.position >= 20]
    target_pam = valid_pams[0] if valid_pams else pams[0]
    
    print(f"🎯 Selected PAM: '{target_pam.sequence}' at position {target_pam.position}")

    try:
        # Simulate cut with NHEJ (creates indel = enhancer disruption)
        result = cas9.cut(
            target=bcl11a_enhancer,
            pam=target_pam,
            repair_pathway="NHEJ",
            cell_type="iPSC"  # Patient-derived stem cells
        )
    except ValueError as e:
        print(f"\n❌ Simulation Failed: {e}")
        return

    print(f"\n--- 📊 RESULTS ---")
    print(f"Cut Efficiency: {result.efficiency:.1f}%")
    print(f"Repair Pathway: {result.repair_pathway}")
    print(f"Cut Position:   {result.cut_position}")

    print(f"\n🧬 Predicted Indel Profile:")
    # Sort indels by probability for better readability
    sorted_indels = sorted(result.indels.items(), key=lambda x: x[1], reverse=True)
    for indel, prob in sorted_indels:
        print(f"  ▸ {indel.replace('_', ' ').title()}: {prob}%")

    print(f"\n⚠️ Off-targets detected: {len(result.off_targets)}")
    for ot in result.off_targets:
        # Updated key from 'cfd_score' to 'simulated_cfd_score' matching previous fixes
        print(f"  - {ot['chromosome']}:{ot['position']} "
              f"(mismatches: {ot['mismatches']}, Simulated CFD: {ot['simulated_cfd_score']})")

    print(f"\n🏥 Clinical Outcome Prediction:")
    # Adjusted thresholds slightly to be realistic with NHEJ profiling
    if result.efficiency > 70 and len(result.off_targets) <= 2:
        print("  ✅ High efficacy predicted. HbF reactivation likely.")
        print("  ✅ Potential functional cure for sickle cell disease.")
    else:
        print("  ⚠️ Optimization needed. Consider alternative gRNA to reduce off-targets or increase efficiency.")
    
    print("=" * 50)

if __name__ == "__main__":
    simulate_sickle_cell_cure()
