@echo off
REM Auto-monitor campaign until completion

echo Starting auto-monitor for campaign...
echo.
echo Campaign is running in background
echo This will monitor until physics validation completes (30 structures)
echo.
echo Press Ctrl+C to stop monitoring (campaign continues in background)
echo.

py -3.12 monitor_until_complete.py

pause
