import { createFileRoute } from "@tanstack/react-router"
import DiazoDB from "../components/DiazoDB"

export const Route = createFileRoute("/classify")({
  component: DiazoDB,
})
