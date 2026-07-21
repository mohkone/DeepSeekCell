from reportlab.lib import colors
from reportlab.lib.pagesizes import landscape
from reportlab.lib.units import inch
from reportlab.pdfgen import canvas
from reportlab.platypus import Paragraph
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.enums import TA_CENTER, TA_LEFT


PAGE = landscape((11 * inch, 6.7 * inch))
OUT = "paper/oup-authoring-template/oup-authoring-template/Fig/deepseekcell_workflow.pdf"


def para(c, text, x, y, w, h, size=9, leading=11, align=TA_CENTER, color=colors.HexColor("#102033")):
    style = ParagraphStyle(
        "box",
        fontName="Helvetica",
        fontSize=size,
        leading=leading,
        alignment=align,
        textColor=color,
        spaceAfter=0,
        spaceBefore=0,
    )
    p = Paragraph(text, style)
    p.wrapOn(c, w, h)
    p.drawOn(c, x, y + (h - p.height) / 2)


def box(c, x, y, w, h, title, body, fill, stroke=colors.HexColor("#35506b")):
    c.setFillColor(fill)
    c.setStrokeColor(stroke)
    c.setLineWidth(1.2)
    c.roundRect(x, y, w, h, 8, stroke=1, fill=1)
    para(c, f"<b>{title}</b>", x + 8, y + h - 25, w - 16, 18, size=10, leading=12)
    para(c, body, x + 9, y + 9, w - 18, h - 36, size=8.3, leading=10)


def arrow(c, x1, y1, x2, y2, color=colors.HexColor("#37536d")):
    c.setStrokeColor(color)
    c.setFillColor(color)
    c.setLineWidth(1.6)
    c.line(x1, y1, x2, y2)
    angle = 0 if x2 >= x1 else 180
    if abs(y2 - y1) > abs(x2 - x1):
        angle = 90 if y2 >= y1 else -90
    import math

    rad = math.radians(angle)
    left = rad + math.radians(155)
    right = rad - math.radians(155)
    size = 7
    p1 = (x2, y2)
    p2 = (x2 + size * math.cos(left), y2 + size * math.sin(left))
    p3 = (x2 + size * math.cos(right), y2 + size * math.sin(right))
    c.line(*p1, *p2)
    c.line(*p1, *p3)


def main():
    c = canvas.Canvas(OUT, pagesize=PAGE)
    width, height = PAGE

    c.setFillColor(colors.white)
    c.rect(0, 0, width, height, fill=1, stroke=0)

    c.setFillColor(colors.HexColor("#0b1f33"))
    c.setFont("Helvetica-Bold", 17)
    c.drawCentredString(width / 2, height - 34, "DeepSeekCell evidence-guided selective refinement")
    c.setFont("Helvetica", 10)
    c.setFillColor(colors.HexColor("#49627a"))
    c.drawCentredString(
        width / 2,
        height - 51,
        "One cached first-pass LLM output is audited by deterministic biological evidence before any second-pass inference is spent.",
    )

    top_y = height - 147
    x0 = 28
    w = 118
    h = 70
    gap = 18

    fills = {
        "input": colors.HexColor("#eaf3ff"),
        "llm": colors.HexColor("#eef8f2"),
        "evidence": colors.HexColor("#fff3df"),
        "decision": colors.HexColor("#f5edff"),
        "refine": colors.HexColor("#ffeceb"),
        "output": colors.HexColor("#eaf7f7"),
    }

    box(c, x0, top_y, w, h, "1. Cluster markers", "Top genes per cluster<br/>tissue and species metadata", fills["input"])
    box(c, x0 + (w + gap), top_y, w, h, "2. First-pass LLM", "Structured JSON labels<br/>confidence, candidates<br/>and reasoning", fills["llm"])
    box(
        c,
        x0 + 2 * (w + gap),
        top_y,
        w + 20,
        h,
        "3. Evidence layer",
        "Cell Ontology mapping<br/>marker-profile support<br/>tissue consistency<br/>mixed-profile signal",
        fills["evidence"],
    )
    box(
        c,
        x0 + 3 * (w + gap) + 20,
        top_y,
        w,
        h,
        "4. Conflict test",
        "Does marker evidence<br/>strongly disagree with<br/>the first-pass label?",
        fills["decision"],
    )
    box(
        c,
        x0 + 4 * (w + gap) + 20,
        top_y,
        w,
        h,
        "5. Refine selected k",
        "Second LLM pass only<br/>for conflicted clusters",
        fills["refine"],
    )

    ymid = top_y + h / 2
    arrow(c, x0 + w, ymid, x0 + w + gap - 5, ymid)
    arrow(c, x0 + 2 * w + gap, ymid, x0 + 2 * (w + gap) - 5, ymid)
    arrow(c, x0 + 2 * (w + gap) + w + 20, ymid, x0 + 3 * (w + gap) + 15, ymid)
    arrow(c, x0 + 3 * (w + gap) + 20 + w, ymid, x0 + 4 * (w + gap) + 15, ymid)

    c.setFont("Helvetica-Bold", 9)
    c.setFillColor(colors.HexColor("#7b2d24"))
    c.drawCentredString(x0 + 4 * (w + gap) + 20 - gap / 2, ymid + 13, "YES")

    # No-conflict branch.
    no_x = x0 + 3 * (w + gap) + 20
    no_y = top_y - 92
    box(
        c,
        no_x - 28,
        no_y,
        w + 56,
        55,
        "No conflict",
        "Retain first-pass label<br/>with evidence-adjusted confidence",
        fills["llm"],
    )
    arrow(c, x0 + 3 * (w + gap) + 20 + w / 2, top_y, no_x + w / 2, no_y + 55)
    c.setFont("Helvetica-Bold", 9)
    c.setFillColor(colors.HexColor("#2d6b45"))
    c.drawString(no_x - 24, no_y + 61, "NO")

    out_x = width - 190
    out_y = top_y - 105
    box(
        c,
        out_x,
        out_y,
        145,
        82,
        "Final annotation",
        "Cell type label<br/>Cell Ontology ID<br/>raw and evidence-adjusted confidence<br/>provenance fields",
        fills["output"],
    )
    arrow(c, x0 + 4 * (w + gap) + 20 + w / 2, top_y, out_x + 80, out_y + 82)
    y_line = out_y + 41
    arrow(c, no_x + w + 28, y_line, out_x, y_line)

    # Benchmark control panel.
    panel_x = 45
    panel_y = 62
    panel_w = width - 90
    panel_h = 118
    c.setFillColor(colors.HexColor("#f7f9fb"))
    c.setStrokeColor(colors.HexColor("#a8b7c6"))
    c.setLineWidth(1)
    c.roundRect(panel_x, panel_y, panel_w, panel_h, 10, stroke=1, fill=1)
    c.setFont("Helvetica-Bold", 12)
    c.setFillColor(colors.HexColor("#102033"))
    c.drawString(panel_x + 16, panel_y + panel_h - 24, "Paired ablation benchmark")
    para(
        c,
        "All DeepSeekCell arms reuse the same cached first-pass response hash. "
        "Plain, Evidence and Calibrated keep labels fixed; selector arms spend matched second-pass budgets.",
        panel_x + 16,
        panel_y + panel_h - 54,
        panel_w - 32,
        28,
        size=9,
        leading=11,
        align=TA_LEFT,
    )

    labels = [
        ("Random-k", "#ffffff"),
        ("Confidence-k", "#ffffff"),
        ("NoOntology-k", "#ffffff"),
        ("Evidence-k", "#fff3df"),
        ("FullRefined", "#ffffff"),
    ]
    chip_w = 128
    chip_gap = 16
    start_x = panel_x + 32
    chip_y = panel_y + 18
    for i, (label, fill) in enumerate(labels):
        x = start_x + i * (chip_w + chip_gap)
        c.setFillColor(colors.HexColor(fill))
        c.setStrokeColor(colors.HexColor("#6b7f93"))
        c.roundRect(x, chip_y, chip_w, 28, 7, stroke=1, fill=1)
        para(c, f"<b>{label}</b>", x + 4, chip_y + 4, chip_w - 8, 20, size=9)

    c.setFont("Helvetica", 8)
    c.setFillColor(colors.HexColor("#5f7284"))
    c.drawRightString(width - 34, 24, "DeepSeekCell 0.1.0")

    c.save()


if __name__ == "__main__":
    main()
