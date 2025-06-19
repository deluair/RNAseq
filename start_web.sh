#!/bin/bash

# RNA-seq Analysis Platform Web Application Starter

echo "🚀 Starting RNA-seq Analysis Platform Web Application"
echo "=================================================="

# Check if port 5000 is available, otherwise use 5001
if ! lsof -Pi :5000 -sTCP:LISTEN -t >/dev/null ; then
    PORT=5000
    echo "✅ Port 5000 is available"
else
    PORT=5001
    echo "⚠️  Port 5000 is in use, using port 5001"
fi

# Kill any existing instances
echo "🔄 Cleaning up any existing processes..."
lsof -ti:5000,5001 | xargs kill -9 2>/dev/null || echo "   No processes to clean up"

# Start the application
echo "🌐 Starting web application on port $PORT..."
echo "   Access at: http://localhost:$PORT"
echo "   Press Ctrl+C to stop"
echo "=================================================="

python app.py --port $PORT 