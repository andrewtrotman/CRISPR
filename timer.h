/*
	TIMER.H
	-------
	Copyright (c) 2017 Andrew Trotman
	Released under the 2-clause BSD license (See:https://en.wikipedia.org/wiki/BSD_licenses)
*/
/*!
	@file
	@brief Timer methods in nanoseconds and milliseconds.
	@author Andrew Trotman
	@copyright 2017 Andrew Trotman
*/
#pragma once

#include <stdio.h>
#include <stdint.h>

#include <chrono>
#include <thread>

namespace JASS
	{
	/*
		CLASS TIMER
		-----------
	*/
	/*!
		@brief General purpose cross-platform timer methods for C++11 and later
	*/
	class timer
		{
		private:
			/*
				TYPEDEF TIMER::STOP_WATCH
				-------------------------
			*/
			/*!
				@typedef stop_watch
				@brief An instance of a timer
			*/
			typedef std::chrono::time_point<std::chrono::steady_clock, std::chrono::nanoseconds> stop_watch;

			/*
				CLASS TIMER::DURATION
				----------------------
			*/
			/*!
				@brief A duration as returned by a stop_watch.
			*/
			class duration
				{
				private:
					std::chrono::nanoseconds took;			///< The time delta that this duration represents.

				public:
					/*
						TIMER::DURATION::DURATION()
						---------------------------
					*/
					/*!
						@brief Convert a std::chrono::nanoseconds into a JASS::timer::duration object.
					*/
					duration(std::chrono::nanoseconds value)
						{
						took = value;
						}

					/*
						TIMER::DURATION::MILLISECONDS()
						-------------------------------
					*/
					/*!
						@brief Return the time-slice as number of milliseconds.
						@return time in miliseconds.
					*/
					auto milliseconds(void)
						{
						return std::chrono::duration_cast<std::chrono::milliseconds>(took).count();
						}

					/*
						TIMER::DURATION::MICROSECONDS()
						-------------------------------
					*/
					/*!
						@brief Return the time-slice as number of microseconds.
						@return time in miliseconds.
					*/
					auto microseconds(void)
						{
						return std::chrono::duration_cast<std::chrono::microseconds>(took).count();
						}

					/*
						TIMER::DURATION::NANOSECONDS()
						------------------------------
					*/
					/*!
						@brief Return the time-slice as number of nanoseconds.
						@return time in nanoseconds.
					*/
					auto nanoseconds(void)
						{
						return std::chrono::duration_cast<std::chrono::nanoseconds>(took).count();
						}
				};

		public:
			/*
				TIMER::START()
				--------------
			*/
			/*!
				@brief Start a stop watfh
				@return A stop watch object that is later passed to timer::stop()
			*/
			static stop_watch start()
				{
				return std::chrono::steady_clock::now();
				}

			/*
				TIMER::STOP()
				-------------
			*/
			/*!
				@brief Return the time on the stop watch.
				@param watch [in] The stop watch to check.
				@return a duration that can be converted into a human-usable unit such as nanoseconds.
			*/
			static duration stop(stop_watch watch)
				{
				return std::chrono::steady_clock::now() - watch;
				}
		};
	}
