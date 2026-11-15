#!/usr/bin/env ruby
# Run an MPI job with a deadline without depending on GNU timeout (e.g. on macOS).
def monotonic_time
  Process.clock_gettime(Process::CLOCK_MONOTONIC)
end

def signal_group(signal, pid)
  Process.kill(signal, -pid)
rescue Errno::ESRCH
end

def terminate_group(pid)
  begin
    signal_group('TERM', pid)
    # Signal handlers only record cancellation, so another signal cannot skip KILL.
    sleep 0.5
  ensure
    signal_group('KILL', pid)
  end
  # Reap without an unbounded wait, even if an OS-level failure delays termination.
  deadline = monotonic_time + 2
  loop do
    begin
      return if Process.waitpid2(pid, Process::WNOHANG)
    rescue Errno::ECHILD
      return
    end
    break if monotonic_time >= deadline
    sleep 0.01
  end
  warn 'Command did not exit after KILL'
end

begin
  seconds = Float(ARGV.shift)
  raise ArgumentError unless seconds.finite? && seconds.positive? && !ARGV.empty?
rescue ArgumentError, TypeError
  warn 'Usage: run_with_timeout.rb POSITIVE_SECONDS COMMAND [ARG ...]'
  exit 2
end

# Record signals rather than raising asynchronously during spawn, cleanup or reap.
cancellation = nil
%w[HUP INT QUIT TERM].each do |signal|
  Signal.trap(signal) { cancellation ||= Signal.list.fetch(signal) }
end
begin
  command = ARGV.shift
  pid = Process.spawn([command, command], *ARGV, pgroup: true)
rescue Errno::ENOENT => error
  warn error.message
  exit 127
end

deadline = monotonic_time + seconds
loop do
  status = Process.waitpid2(pid, Process::WNOHANG)&.last
  if cancellation
    terminate_group(pid)
    exit 128 + cancellation
  end
  exit(status.exitstatus || 128 + status.termsig) if status
  if monotonic_time >= deadline
    terminate_group(pid)
    warn "Command timed out after #{seconds} seconds"
    exit 124
  end
  sleep 0.01
end
