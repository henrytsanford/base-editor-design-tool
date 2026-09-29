"""The base-editor guide designer, exposed over the Model Context Protocol.

Imports the engine directly, the same way service/jobs.py does, and reads the same
local reference bundle. What the model gets is the web app's surface without the
HTML: gene and transcript lookup, a design run, and filtered views of its result.

Run it with `python -m mcp_server`. Nothing is imported here, so that module can be
executed without being imported twice.
"""
