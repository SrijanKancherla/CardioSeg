import numpy as np 
from cellpose import models
from skimage.measure import label
from skimage.morphology import remove_small_objects

class CellposeSegmenter:
    def __init__(self, gpu=True):
        from cellpose import models
        self.model = models.CellposeModel(
            gpu=gpu
            )

    def segment(
        self,
        image,
        diameter=None,
        progress_callback=None,
        flow_threshold=0.4,
        cellprob_threshold=0.0,
        min_size=15,
    ):
        """
        Segment nuclei by unioning:
        - high-recall pass: flow_threshold=0.4, cellprob_threshold=0.0
        - high-precision pass: flow_threshold=0.8, cellprob_threshold=-1.0
        with min_size=5 in both passes.
        """
        recall_masks, _, _ = self._run_eval(
            image=image,
            diameter=diameter,
            flow_threshold=0.4,
            cellprob_threshold=0.0,
            min_size=5,
            progress_callback=progress_callback,
            progress_start=0.0,
            progress_end=0.5,
        )


        precision_masks, flows, styles = self._run_eval(
            image=image,
            diameter=diameter,
            flow_threshold=0.8,
            cellprob_threshold=-1.0,
            min_size=5,
            progress_callback=progress_callback,
            progress_start=0.5,
            progress_end=1.0,
        )

        masks = self._instance_merge(recall_masks, precision_masks)

        if progress_callback is not None:
            progress_callback(1.0)

        return {
            "masks": masks,
            "flows": flows,
            "cellprob": flows[2] if len(flows) > 2 else None
        }
    
    def _run_eval(
        self,
        image,
        diameter,
        flow_threshold,
        cellprob_threshold,
        min_size,
        progress_callback,
        progress_start,
        progress_end,
    ):
        scaled_progress = None
        if progress_callback is not None:
            span = progress_end - progress_start

            def scaled_progress(raw_value):
                progress_callback(progress_start + span * float(raw_value))

        masks, flows, styles = self.model.eval(
            image,
            diameter=diameter,
            progress=scaled_progress,
            augment=False,
            normalize=True,
            flow_threshold=flow_threshold,
            cellprob_threshold=cellprob_threshold,
            min_size=min_size,
        )
        return masks, flows, styles

    def _instance_merge(
        self,
        mask_a,
        mask_b,
        iou_thresh=0.5,
        containment_thresh=0.85,
        min_size=5,
    ):
        """
        Improved instance merge:
        - Matches precision → recall (IoU)
        - Removes nested masks (containment)
        - Keeps unmatched recall masks (rescues)
        """

        def get_instances(mask):
            ids = np.unique(mask)
            ids = ids[ids != 0]
            return {i: (mask == i) for i in ids}

        inst_a = get_instances(mask_a)  # recall
        inst_b = get_instances(mask_b)  # precision

        used_a = set()
        merged = []

        # --- Step 1: match precision → recall ---
        for id_b, m_b in inst_b.items():
            best_iou = 0
            best_a = None
            best_containment = 0

            for id_a, m_a in inst_a.items():
                intersection = np.logical_and(m_a, m_b).sum()
                if intersection == 0:
                    continue

                union = np.logical_or(m_a, m_b).sum()
                iou = intersection / union

                containment = intersection / m_b.sum()  # how much of B is inside A

                if iou > best_iou:
                    best_iou = iou
                    best_a = id_a
                    best_containment = containment

            # --- Decision logic ---
            if best_a is not None and (
                best_iou > iou_thresh or best_containment > containment_thresh
            ):
                # Merge into recall object
                merged_mask = np.logical_or(inst_a[best_a], m_b)
                used_a.add(best_a)
                merged.append(merged_mask)
            else:
                # Keep precision object
                merged.append(m_b)

        # --- Step 2: add unmatched recall (rescued cells) ---
        for id_a, m_a in inst_a.items():
            if id_a not in used_a:
                merged.append(m_a)

        # --- Step 3: remove nested masks ---
        def remove_nested(masks):
            keep = []
            for i, m1 in enumerate(masks):
                area1 = m1.sum()
                is_nested = False

                for j, m2 in enumerate(masks):
                    if i == j:
                        continue

                    intersection = np.logical_and(m1, m2).sum()
                    if intersection == 0:
                        continue

                    containment = intersection / area1

                    if containment > containment_thresh:
                        is_nested = True
                        break

                if not is_nested:
                    keep.append(m1)

            return keep

        merged = remove_nested(merged)

        # --- Step 4: rebuild label image ---
        final = np.zeros(mask_a.shape, dtype=np.uint32)
        label_id = 1

        for m in merged:
            if m.sum() >= min_size:
                final[m] = label_id
                label_id += 1

        return final


#What is flows?
'''
flows[0] -> dY flow, flows[1] -> dX flow, flows[2] -> cell probability map, flows[3] -> internal pixel
'''

