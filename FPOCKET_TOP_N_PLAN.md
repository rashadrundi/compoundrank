# fpocket top-N implementation plan

Goal:
Allow automatic pocket mode to dock the top N fpocket pockets instead of only pocket 1.

Why:
Single-pocket fpocket worked technically but selected a non-catalytic HIV protease pocket. Top-N pocket docking reduces the chance of missing the true active site.

Rules:
- Explicit box mode stays unchanged.
- Autobox ligand mode stays unchanged.
- fpocket mode becomes exploratory multi-pocket mode.
- No custom pose ranking is added.
- GNINA score remains the only ordering source.
- Final PDBs must label which pocket was used.
- Temporary fpocket paths must not appear in final REMARK 900 output.

Needed changes:
1. Add CLI option:
   --fpocket-top-n

2. Update pocket.py:
   fpocket should return one or more PocketDefinition objects.

3. Update pipeline.py:
   if multiple pockets exist, dock each ligand against each pocket.

4. Update export.py / output naming:
   include pocket identity in REMARK 900 and possibly filename.

5. Add tests:
   verify pocket definitions can carry pocket labels.
   verify multi-pocket output naming does not collide.
