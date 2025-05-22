import webbrowser
import os

class BioinfoCalculator:
    def __init__(self):
        self.total_score = 0
        self.total_credits = 0

    def input_scores(self):
        print("🎓 Welcome to the MSc Bioinformatics Grade Calculator!\n")

        # Dissertation status
        self.completed_diss = input("📘 Have you completed your dissertation? (yes/no): ").strip().lower()

        # Optional unit type
        self.optional_unit = input(
            "📚 Which optional unit did you take?\n   👉 Applied Data Science in Biology\n   👉 Molecular Phylogenetics and Epidemiology\nYour choice: "
        ).strip().lower()

        # --- Advances --- #
        print("\n🧠 Enter your marks for *Advances in Bioinformatics* assignments (as % e.g. 75):")
        a1 = float(input("   📝 Assignment 1 (20%): "))
        a2 = float(input("   📝 Assignment 2 (20%): "))
        a3 = float(input("   📝 Assignment 3 (30%): "))
        a4 = float(input("   📝 Assignment 4 (30%): "))
        self.advances_score = a1 * 0.2 + a2 * 0.2 + a3 * 0.3 + a4 * 0.3

        # --- Optional Unit ---
        print(f"\n📦 Enter your marks for the *{self.optional_unit.upper()}* optional unit:")
        if "data" in self.optional_unit:
            c1 = float(input("   📊 Coursework 1 (40%): "))
            c2 = float(input("   📄 Report (60%): "))
            self.optional_score = c1 * 0.4 + c2 * 0.6
        elif "molecular" in self.optional_unit:
            m1 = float(input("   🧬 Assessment 1 (20%): "))
            m2 = float(input("   🧬 Assessment 2 (20%): "))
            m3 = float(input("   🧪 Exam (60%): "))
            self.optional_score = m1 * 0.2 + m2 * 0.2 + m3 * 0.6
        else:
            print("❌ Invalid optional unit selected. Please restart and check your input.")
            return False

        # --- RP1b --- #
        print("\n🔬 Enter your marks for *Research Project 1B*:")
        r1 = float(input("   🧪 Assessment 1 (20%): "))
        r2 = float(input("   🧪 Assessment 2 (20%): "))
        r3 = float(input("   📄 Report (60%): "))
        self.rp1b_score = r1 * 0.2 + r2 * 0.2 + r3 * 0.6

        # --- BHOR --- #
        print("\n🌍 Enter your marks for *Broadening Horizons*:")
        pres = float(input("   🎤 Presentation/Engagement (50%): "))
        peer = float(input("   🤝 Peer Mark (15%): "))
        portfolio = float(input("   📁 Individual Portfolio (35%): "))
        self.bhor_score = pres * 0.5 + peer * 0.15 + portfolio * 0.35

        # --- RP2 --- #
        if self.completed_diss == "yes":
            self.rp2_score = float(input("\n📘 Enter your *Dissertation (RP2)* mark (100%): "))
        else:
            self.rp2_score = None

        return True

    def calculate_average(self):
        advances_weighted = self.advances_score * 40
        optional_weighted = self.optional_score * 30
        rp1b_weighted = self.rp1b_score * 30
        bhor_weighted = self.bhor_score * 20

        self.total_score = advances_weighted + optional_weighted + rp1b_weighted + bhor_weighted

        if self.rp2_score is not None:
            rp2_weighted = self.rp2_score * 60
            self.total_score += rp2_weighted
            self.total_credits = 180
        else:
            self.total_credits = 120

        average = self.total_score / self.total_credits
        print(f"\n✅ Your weighted Bioinformatics average is: 🎯 {average:.2f}%")
        return average

    def generate_html_report(self, average):
        html_content = f"""
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>MSc Bioinformatics Grade Report</title>
            <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;600&display=swap" rel="stylesheet">
            <style>
                body {{
                    font-family: 'Inter', sans-serif;
                    background: #f0f4f8;
                    margin: 0;
                    padding: 2rem;
                    color: #333;
                }}
                .container {{
                    max-width: 700px;
                    margin: 0 auto;
                    background: white;
                    border-radius: 12px;
                    box-shadow: 0 4px 20px rgba(0, 0, 0, 0.08);
                    padding: 2rem 2.5rem;
                }}
                h1 {{
                    text-align: center;
                    color: #2E8B57;
                    margin-bottom: 1.5rem;
                }}
                .section {{
                    margin-bottom: 1.5rem;
                }}
                .label {{
                    font-weight: 600;
                    color: #555;
                }}
                .value {{
                    float: right;
                    font-weight: 600;
                    color: #000;
                }}
                .score {{
                    display: flex;
                    justify-content: space-between;
                    padding: 0.5rem 0;
                    border-bottom: 1px solid #eee;
                }}
                .final {{
                    font-size: 1.4rem;
                    text-align: center;
                    color: #2E8B57;
                    font-weight: 700;
                    margin-top: 2rem;
                    padding-top: 1rem;
                    border-top: 2px solid #eee;
                }}
            </style>
        </head>
        <body>
            <div class="container">
                <h1>📘 MSc Bioinformatics Grade Report</h1>

                <div class="section">
                    <div class="score"><span class="label">Dissertation completed:</span><span class="value">{self.completed_diss.capitalize()}</span></div>
                    <div class="score"><span class="label">Optional Unit:</span><span class="value">{self.optional_unit.title()}</span></div>
                </div>

                <div class="section">
                    <div class="score"><span class="label">Advances Score (x40):</span><span class="value">{self.advances_score:.2f}</span></div>
                    <div class="score"><span class="label">Optional Unit Score (x30):</span><span class="value">{self.optional_score:.2f}</span></div>
                    <div class="score"><span class="label">RP1b Score (x30):</span><span class="value">{self.rp1b_score:.2f}</span></div>
                    <div class="score"><span class="label">BHOR Score (x20):</span><span class="value">{self.bhor_score:.2f}</span></div>"""

        if self.rp2_score is not None:
            html_content += f"""
                    <div class="score"><span class="label">RP2 Dissertation Score (x60):</span><span class="value">{self.rp2_score:.2f}</span></div>"""

        html_content += f"""
                </div>

                <div class="final">
                    🎯 Final Weighted Average: {average:.2f}%
                </div>
            </div>
        </body>
        </html>
        """

        filename = "bioinfo_report.html"
        with open(filename, "w") as f:
            f.write(html_content)

        print(f"📄 HTML report generated: {filename}")
        webbrowser.open(f"file://{os.path.abspath(filename)}")


# --- 🎬 Example Usage --- #
if __name__ == "__main__":
    calculator = BioinfoCalculator()
    if calculator.input_scores():
        avg = calculator.calculate_average()
        calculator.generate_html_report(avg)
