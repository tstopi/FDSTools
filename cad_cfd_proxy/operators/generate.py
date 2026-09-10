# -*- coding: utf-8 -*-
"""Generate operator — runs the pipeline with progress and cancel.

Interactive invoke runs a modal, timer-driven loop that advances the pipeline
one phase per tick, updates a progress bar, and honours cancel between phases
(the heavy per-phase ops can't be interrupted mid-step). Non-interactive
execute (scripts, ``bpy.ops`` without a window) falls back to a synchronous run.
"""

import bpy
from bpy.types import Operator

from .. import compat, core

_TIMER_STEP = 0.05  # seconds between pipeline phases


class CADCFD_OT_generate(Operator):
    bl_idname = "cadcfd.generate"
    bl_label = "Generate Envelope"
    bl_description = "Build a simplified CFD proxy mesh from the source assembly"
    bl_options = {"REGISTER", "UNDO"}

    _timer = None
    _gen = None

    @classmethod
    def poll(cls, context):
        missing = compat.missing_capabilities()
        if missing:
            cls.poll_message_set("Unsupported Blender build: " + ", ".join(missing))
            return False
        props = context.scene.cad_cfd_proxy
        if props.is_running:
            cls.poll_message_set("A generation is already running")
            return False
        return True

    # --- interactive (modal, with progress + cancel) ---

    def invoke(self, context, event):
        props = context.scene.cad_cfd_proxy
        # Creating the generator does not run any phase (lazy); the first phase
        # runs on the first modal tick, where errors are handled.
        self._gen = core.generate_job(context, props)

        props.is_running = True
        props.cancel_requested = False
        props.progress = 0.0
        props.progress_label = "Starting"

        wm = context.window_manager
        wm.progress_begin(0.0, 1.0)
        self._timer = wm.event_timer_add(_TIMER_STEP, window=context.window)
        wm.modal_handler_add(self)
        return {"RUNNING_MODAL"}

    def modal(self, context, event):
        props = context.scene.cad_cfd_proxy

        if event.type == "ESC" or props.cancel_requested:
            return self._end(context, cancelled=True)

        if event.type != "TIMER":
            return {"RUNNING_MODAL"}

        try:
            label, fraction = next(self._gen)
        except StopIteration:
            return self._end(context, cancelled=False)
        except Exception as exc:  # noqa: BLE001
            # Any phase error (expected PipelineError, or an unexpected C-op
            # failure) must still tear the modal down — otherwise is_running
            # stays True and the event timer leaks, wedging Generate.
            self.report({"ERROR"}, str(exc))
            return self._end(context, cancelled=True)

        props.progress = fraction
        props.progress_label = label
        context.window_manager.progress_update(fraction)
        _tag_redraw(context)
        return {"RUNNING_MODAL"}

    def _end(self, context, cancelled):
        props = context.scene.cad_cfd_proxy
        wm = context.window_manager
        if self._timer is not None:
            wm.event_timer_remove(self._timer)
            self._timer = None
        wm.progress_end()
        self._gen = None
        props.is_running = False
        props.cancel_requested = False
        props.progress = 0.0
        props.progress_label = ""
        _tag_redraw(context)
        if cancelled:
            core.discard_work(context)
            self.report({"WARNING"}, "Generation cancelled")
            return {"CANCELLED"}
        self.report({"INFO"}, "Proxy generated")
        return {"FINISHED"}

    # --- non-interactive fallback (scripts / headless) ---

    def execute(self, context):
        props = context.scene.cad_cfd_proxy
        try:
            core.generate_proxy(context, props, report=self.report)
        except (NotImplementedError, core.PipelineError) as exc:
            core.discard_work(context)  # don't leave partial intermediates behind
            self.report({"ERROR"}, str(exc))
            return {"CANCELLED"}
        self.report({"INFO"}, "Proxy generated")
        return {"FINISHED"}


def _tag_redraw(context):
    if context.area is not None:
        context.area.tag_redraw()
